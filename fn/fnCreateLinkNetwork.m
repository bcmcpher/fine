function [ netw ] = fnCreateLinkNetwork(netw, meas, dtype, Phi, dict, clobber)
%fnCreateLinkNetwork creates a link network from a pconn list that has had
% edge volumes precomputed for the label requested.
%
% INPUTS:
%     pconn - is the paired connections object to create a link network from. 
%             Edge volumes need to be estimated beforehand with fnFindPathVoxels.m 
%
%     meas  - the measure to be used when creating the links. 
%             Possible options are:
%               - dice: the dice coefficient between voxels of the connections
%               - mi: mutual information between precomputed voxel values
%               - angle: the average angle of intersection between the streamline nodes
%               - tprof: a series of summary statistics between a specified tract profile
%
%     dtype - defines the microstructural label for "mi" or the profile label for "tprof"
%
%     Phi   - the Phi model object from the fe structure. 
%             Only needed if "angle" requested.
%
%     dict  - the dictionary of streamline orientations from the fe structure. 
%             Only needed if "angle" requested.
%
% OUTPUTS:
%     omat - a matrix that is (length(pconn) x length(pconn)) in size
%            showing the requested measure of each intersecting edge.
%
%     olab - a cell array of labels that describe the 3rd dimension of omat
%
%     out  - cell array of the data used to compute edge entry in omat
%
% TODO:
% - add other metrics?
%   https://en.wikipedia.org/wiki/Diversity_index#Simpson_index
%     - Simpson's Index: probability that 2 random voxels are part of intersection
%     - Richness: how many voxels does an intersection contain
%     - Shannon Diversity Index: global measure of intersections between links 
%
% EXAMPLE:
%
% % load data
% parc          = niftiRead('labels.nii.gz');
% fg            = feGet(fe, 'fibers acpc');
% fibers        = fg.fibers;
% fibLength     = fefgGet(fg, 'length');
% weights       = feGet(fe, 'fiberweights');
% Phi           = feGet(fe, 'Phi');
% vx            = [ 1.25 1.25 1.25 ];
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % find the voxels that each edge occupies
% pconn = fnFindPathVoxels(pconn, 'nzw', Phi, vx);
%
% % create link network based on overlap of non-zero weighted fibers
% [ omat, out ] = fnCreateLinkNetwork(pconn, 'nzw', 'dice');
%
% Brent McPherson (c), 2017 - Indiana University
%

%% generate link network

% if dtype isn't passed, set as empty b/c it isn't used
if(~exist('dtype', 'var') || isempty(dtype))
    dtype = 'none';
end

% if clobber isn't passed, set as false
if(~exist('clobber', 'var') || isempty(clobber))
    clobber = false;
end

% if a link network is already computed and clobber is true, throw warning
if isfield(netw, 'links') && clobber
    warning('Replacing previously computed link network with requested values.');
end

% if a link network is already computed and clobber is false, throw error
if isfield(netw, 'links') && ~clobber
    error('Link network is already computed and clobber is set to false.'); 
end

% error if coordinates are not stored from volume averaging
if ~isfield(netw.edges{1}.volume, 'coords')
    error('Edge volume coordinates must be precomputed and stored with fnAveragePropertyEdges()');
end

% error if MI is called to estimate w/o dtype being stored
if strcmp(meas, 'mi') && ~isfield(netw.edges{1}, 'micros')
    if ~isfield(netw.edges{1}.micros, dtype)
        if ~isfield(netw.edges{1}.micros.(dtype), 'raw')
            error('An average edge property dtype ''%s'' must be precomputed with fnAveragePropertyEdges()', dtype);
        end
    end
end

% error if tract profiles are called w/o dtype being stored
if strcmp(meas, 'tprof') && ~isfield(netw.edges{1}, 'profile')
    if ~isfield(netw.edges{1}.profile, dtype)
        error('Tract profiles for ''%s'' must be precomputed and stored with fnTractProfileEdges().', dtype);
    end
end

% error when Phi needed and not found
if(~exist('Phi', 'var') || isempty(Phi))
    Phi = [];
    if strcmp(meas, 'angle')
        error('Phi was not passed, so the angle cannot be computed.');
    end 
end

% error when dict needed and not found
if(~exist('dict', 'var') || isempty(dict))
    dict = [];
    if strcmp(meas, 'angle')
        error('dict was not passed, so the angle cannot be computed.');
    end 
end

% parese requested measurement to preallocate output
switch meas    
    case 'dice'
        olab = {'Dice Coefficient'};
        prnt = 'dice coefficient';
        outd = 3;
    case 'mi'
        olab = {'Mutual Information', 'Joint Entropy', 'Entropy1', 'Entropy2'};
        nbin = 256; % number of bins for the individual / joint histograms
        prnt = [ 'mutual information of ' dtype ];
        outd = 6;
    case 'angle'
        olab = {'Angle of Intersection', 'Product', 'Cosine Distance'};
        prnt = 'angle of intersection';
        outd = 5;
    case 'tprof'
        olab = {'Norm', 'Correlation', 'Sum of Absolute Error', 'Product', 'RMSE', 'Cosine Distance'};
        prnt = [ 'tract profile similarities of ' dtype ];
        outd = 7;
    otherwise
        error('Invalid method of comparison requested. Either use: ''dice'', ''mi'', ''angle'', ''tprof'', or ''deReus''.');
end

disp([ 'Computing ' prnt ' for unique edge intersections.' ]);

disp('Determining output size...');

% pull the number of nodes to rebuild possible edge list
nnodes = size(netw.nodes, 1);

% build the upper diagonal node indices
nmsk = triu(ones(nnodes), 1);
[ xnidx, ynidx ] = find(nmsk > 0);

% build the possible edges that could exist
eindx = sortrows([ xnidx, ynidx ]);

clear nnodes nmsk xnidx ynidx

% pull the observed (real) edge indices
pairs = netw.parc.pairs;

% pull the possible comparisons to make from edges that are stored
plnk = triu(ones(size(pairs, 1)), 1);
[ xpidx, ypidx ] = find(plnk > 0);

% build the indices into the edges so only all real comparisons are made
tindx = sortrows([ xpidx, ypidx ]);

clear plnk xpidx ypidx

% pull the count of total links to test
nlinks = size(tindx, 1);

% build an empty output array of the largest possible size
out = nan(nlinks, outd);

% pull edges to index the values stored internally
edges = netw.edges;

disp([ 'Testing a possible ' num2str(nlinks) ' edge pairs for intersection...' ]);

tic;
for link = 1:nlinks
    
    % grab the edge indices for the link
    eidx1 = tindx(link, 1);
    eidx2 = tindx(link, 2);

    % grab the unique voxels of each edge
    pv1 = edges{eidx1}.volume.coords;
    pv2 = edges{eidx2}.volume.coords;

    % find the common volume of the connections
    vind = intersect(pv1, pv2, 'rows');
    
    % if the intersection is empty, fill in zero
    if isempty(vind)
        continue
    end
    
    % recover edge index for real edges from all possible edges
    [ ~, edge1 ] = intersect(eindx, pairs(eidx1, :), 'rows');
    [ ~, edge2 ] = intersect(eindx, pairs(eidx2, :), 'rows');
    
    switch meas
        
        % dice coefficient
        case 'dice'
        
            % find the size of the intersection - numerator
            num = size(vind, 1);
            
            % combine the voxel indices for denominator of Dice coeff
            den = size(pv1, 1) + size(pv2, 1);
            
            % build Dice coeff values for assignment
            out(link, :) = [ edge1, edge2, ((2 * num) / den) ];
        
        case 'mi'

            % pull the stored values
            imv1 = edges{eidx1}.micros.(dtype).raw;
            imv2 = edges{eidx2}.micros.(dtype).raw;
            
            % grab combined voxels from both paths
            imp = union(pv1, pv2, 'rows'); 
            
            % build image vector
            im = nan(size(imp, 1), 2);
            
            % grab indices of voxels not in the path
            [ ~, im01 ] = setdiff(imp, pv1, 'rows');
            [ ~, im02 ] = setdiff(imp, pv2, 'rows');
            
            % set entries unique to each tract to zero in the other tract
            im(im01, 1) = 0;
            im(im02, 2) = 0;
            
            % replace missing in each column with the values
            im(isnan(im(:, 1)), 1) = imv1;
            im(isnan(im(:, 2)), 2) = imv2;

            % create first edge histogram and compute entropy
            ihist1 = histcounts(im(:, 1), nbin);
            ihist1 = ihist1 / sum(ihist1);
            ihist1nz = ihist1(ihist1 ~= 0);
            entropy1 = -sum(ihist1nz(:) .* log2(ihist1nz(:)));
            
            % create second edge histogram and compute entropy
            ihist2 = histcounts(im(:, 2), nbin);
            ihist2 = ihist2 / sum(ihist2);
            ihist2nz = ihist2(ihist2 ~= 0);
            entropy2 = -sum(ihist2nz(:) .* log2(ihist2nz(:)));
            
            % create joint edge histogram and compute joint entropy
            jhist = histcounts2(im(:, 1), im(:, 2), nbin);
            jhist = jhist / sum(jhist(:));
            jhistnz = jhist(jhist ~=0);
            jointEntropy = -sum(jhistnz(:) .* log2(jhistnz(:)));
            
            % compute mutual information
            mutualInfo = (entropy1 + entropy2) - jointEntropy;

            % assign to output
            out(link, :) = [ edge1, edge2, mutualInfo, jointEntropy, entropy1, entropy2 ];
            
        % angle of intersection (subset the tensor)    
        case 'angle'
            
            % grab the streamline indices of each edge
            pi1 = edges{eidx1}.fibers.indices;
            pi2 = edges{eidx2}.fibers.indices;
            
            % reduce Phi to the shared streamlines of the connections
            sub1 = Phi(:, :, pi1);
            sub2 = Phi(:, :, pi2);
            
            % from subtensors pull the dictionary atoms for each edge
            atom1 = dict(:, sub1.subs(:, 1));
            atom2 = dict(:, sub2.subs(:, 1));
            
            % find the average orientation of each tract
            matm1 = mean(atom1, 2);
            matm2 = mean(atom2, 2);
            
            % find the scalar product to assign a link weight
            dotp = dot(matm1, matm2);
            
            % compute the angle between connections
            angl = rad2deg(acos(dotp));
            
            % compute cosine distance between angle
            cosd = pdist2(matm1', matm2', 'cosine');
            
            out(link, :) = [ edge1, edge2, angl, dotp, cosd ];
            
        case 'tprof'
            
            % pull tract profiles
            tp1 = edges{eidx1}.profile.(dtype);
            tp2 = edges{eidx2}.profile.(dtype);
            
            % compute the norm between profiles
            nrm = norm(tp1 - tp2);
            
            % other measures
            corc = corr(tp1, tp2);               % correlation
            smae = norm(tp1 - tp2, 1);           % sum of absolute error
            dotp = dot(tp1, tp2);                % scalar product
            rmse = sqrt(mean((tp1 - tp2) .^ 2)); % RMSE
            cosd = pdist2(tp1', tp2', 'cosine'); % cosine distance
            
            out(link, :) = [ edge1, edge2, nrm, corc, smae, dotp, rmse, cosd ];

        otherwise
            error('Invalid metric between edges requested.');            
    end
    
end
time = toc;

disp(['Computed link network metrics in ' num2str(time/60) ' minutes.' ]);

%% store in network object

% drop the empty observations
out(isnan(out(:,3)),:) = [];

disp([ 'Storing ' num2str(size(out, 1)) ' unique links between edges.' ]);

% only store count, non-zero values and labels
netw.links.meas = meas;
netw.links.dtype = dtype;
netw.links.labels = [ 'Edge1', 'Edge2', olab ];
netw.links.values = out;

end
