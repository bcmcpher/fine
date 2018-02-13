function [ omat, olab, out ] = devLinkNetwork(pconn, label, meas, dtype, Phi, dict)
%devLinkNetwork Summary of this function goes here
%   Detailed explanation goes here
%
% TODO:
% - documentation
% - argument parsing / error checking of inputs
%
%

%% generate link network

% error if volume is not precomputed
if (~isfield(pconn{1}.(label), 'volume'))
    error('Edge volumes for this label must be precomputed with fnFindPathVoxels().');
end    

% error if MI is called to estimate w/o dtype being stored
if (~isfield(pconn{1}.(label).volume, 'dtype') && strcmp(meas, 'mi'))
    error('An average edge property for this dtype must be precomputed with fnAverageEdgeProperty().');
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

% check if dtype exists and can be split on a '.', meaning it's a dt6.(type) profile
if(~exist('dtype', 'var') || isempty(dtype))
    dtype = 'none';
end
dtype = strsplit(dtype, '.');

% make this a switch if dtype is needed

disp('Preallocating output values...');

% number of edges; links between edges are the new nodes
lnodes = size(pconn, 1);

% compute indices of output size
mask = triu(ones(lnodes), 1); 

% preallocate unique indices - faster than nchoosek
[ xind, yind ] = find(mask > 0);
indx = [ xind, yind ];

clear mask xind yind 

switch meas
    
    case 'dice'
        olab = {'dice'};
        nmiss = 1;
    case 'mi'
        olab = {'mi', 'joint', 'entropy1', 'entropy2'};
        nmiss = 4;
    case 'angle'
        olab = {'dotp', 'angle', 'cosd'};
        nmiss = 3;
    case 'tprof'
        olab = {'nrm', 'corr', 'smae', 'dotp', 'rmse', 'cosd'};
        nmiss = 6;
    case 'deReus'
        olab = {'node'};
        nmiss = 1;
    otherwise
        error('Invalid method of comparison requested. Either use: ''dice'', ''mi'', ''angle'', ''tprof'', or ''deReus''.');
end

% build an empty output array
out = cell(size(indx, 1), 1);

disp('Computing weights for unique edge intersections...');

tic;
parfor ii = 1:size(indx, 1)
    
    % grab the precomputed indices
    conn1 = indx(ii, 1);
    conn2 = indx(ii, 2);

    % grap the unique streamlines of each edge
    pi1 = pconn{conn1}.(label).indices;
    pi2 = pconn{conn2}.(label).indices;
    
    % fill in nothing if either connection is empty
    if isempty(pi1) || isempty(pi2)
        out{ii} = zeros(1, nmiss);
        continue
    end
    
    % grab the unique voxels of each edge
    pv1 = pconn{conn1}.(label).volume.pvoxels;
    pv2 = pconn{conn2}.(label).volume.pvoxels;
        
    % find the common volume of the connections
    vind = intersect(pv1, pv2);
    
    % if the intersection is empty, fill in zero
    if isempty(vind)
        out{ii} = zeros(1, nmiss);
        continue
    end
    
    switch meas
        
        % dice coefficient
        case 'dice'
        
            % find the size of the intersection - numerator
            num = size(vind, 1);
            
            % combine the voxel indices for denominator of Dice coeff
            den = size(pv1, 1) + size(pv2, 1);
            
            % build Dice coeff values for assignment
            out{ii} = (2 * num) / den;
        
        case 'mi'
            
            % SOME INDIVIDUAL ENTROPY ESTIMATES ARE LARGE NEGATIVE NUMBERS
            % THIS MAKES ~HALF OF MI ESTIMATES LARGE NEGATIVE NUMBERS - BAD
            % I DO NOT CURRENTLY KNOW HOW TO FIX THIS
            % link to how I compute the values currently
            % https://stackoverflow.com/questions/23691398/mutual-information-and-joint-entropy-of-two-images-matlab
            
            % grab path voxels and values for the first connection
            imp1 = pconn{conn1}.(label).volume.pvoxels;
            imv1 = pconn{conn1}.(label).volume.(dtype).raw;
            
            % grab path voxels and values for the second connection
            imp2 = pconn{conn2}.(label).volume.pvoxels;
            imv2 = pconn{conn2}.(label).volume.(dtype).raw;
            
            % grab combined voxels from both paths
            imp = union(imp1, imp2);
            
            % build image vector
            im = nan(size(imp, 1), 2);
            
            % grab indices of voxels not in the path
            [ ~, im01 ] = setdiff(imp, imp1);
            [ ~, im02 ] = setdiff(imp, imp2);
            
            % set entries unique to each tract to zero in the other tract
            im(im01, 1) = 0;
            im(im02, 2) = 0;
            
            % replace missing in each column with the values
            im(isnan(im(:, 1)), 1) = imv1;
            im(isnan(im(:, 2)), 2) = imv2;

            % compute bins of joint histogram
            [ ~, ~, indrow ] = unique(im(:, 1));
            [ ~, ~, indcol ] = unique(im(:, 2));
            
            % compute joint entropy
            jointHist = accumarray([indrow indcol], 1);
            jointProb = jointHist / numel(indrow);
            indNoZero = jointHist ~= 0;
            jointNzPb = jointProb(indNoZero);
            jointEntropy = -sum(jointNzPb .* log2(jointNzPb));
            
            % compute individual histogram summaries
            histImage1 = sum(jointHist, 1);
            histImage2 = sum(jointHist, 2);
            
            % find non-zero elements for first image's histogram
            % extract them out and get the probabilities
            % compute the entropy
            indNoZero1 = histImage1 ~= 0;
            prob1NoZero = histImage1(indNoZero1) / numel(histImage1);
            entropy1 = -sum(prob1NoZero .* log2(prob1NoZero));
            
            % repeat for the second image
            indNoZero2 = histImage2 ~= 0;
            prob2NoZero = histImage2(indNoZero2) / numel(histImage2);
            entropy2 = -sum(prob2NoZero .* log2(prob2NoZero));
            
            % now compute mutual information
            mutualInfo = (entropy1 + entropy2) - jointEntropy;
 
            % assign to output
            out{ii} = [ mutualInfo, jointEntropy, entropy1, entropy2 ];
            
        % angle of intersection (subset the tensor)    
        case 'angle'
            
            % reduce Phi to the shared voxels of the connections
            sub1 = Phi(:, vind, pi1);
            sub2 = Phi(:, vind, pi2);
            
            % from subtensors pull the dictionary atoms for each connection
            atom1 = dict(:, sub1.subs(:, 1));
            atom2 = dict(:, sub2.subs(:, 1));
            
            % find the average orientation of each tract
            matm1 = mean(atom1, 2);
            matm2 = mean(atom2, 2);
            
            % find the scalar product to assign a link weight
            dotp = dot(matm1, matm2);
            
            % compute the angle between connections
            angl = acos(dotp);
            
            % compute cosine distance between angle
            cosd = pdist2(matm1', matm2', 'cosine');
            
            out{ii} = [ dotp angl cosd ];
            
        case 'tprof'
            
            if size(dtype, 2) == 1
                tp1 = pconn{conn1}.(label).profile.(dtype{1});
                tp2 = pconn{conn2}.(label).profile.(dtype{1});
            elseif size(dtype, 2) == 2
                tp1 = pconn{conn1}.(label).profile.(dtype{1}).(dtype{2});
                tp2 = pconn{conn2}.(label).profile.(dtype{1}).(dtype{2});
            else
                error('Impossible tract profile requested.');
            end
            
            % compute the norm between profiles
            nrm = norm(tp1 - tp2);
            
            % other measures
            corc = corr(tp1, tp2);               % correlation
            %mnae = mae(tp1 - tp2);               % mean of absolute error
            smae = norm(tp1 - tp2, 1);           % sum of absolute error
            %mxae = norm(tp1 - tp2, inf);         % max of absolute error 
            dotp = dot(tp1, tp2);                % scalar product
            %mser = mse(tp1, tp2);                % mean of square error
            %sser = sse(tp1 - tp2);               % sum of square error
            rmse = sqrt(mean((tp1 - tp2) .^ 2)); % RMSE
            cosd = pdist2(tp1', tp2', 'cosine'); % cosine distance
            
            % combine and fit coefficients to the lines
            % how could these be compared?
            %[ coeff, errs ] = lpc([ tp1, tp2 ], 4);
            
            % fourier? other time series analysis?
            
            % store output
            %out{ii} = [ nrm corc mnae smae mxae scal mser sser rmse cosd ];
            out{ii} = [ nrm corc smae dotp rmse cosd ];
            
        case 'deReus'
            
            % assume empty
            val = 0;
            
            % if two connections share a node, assign 1
            if (pconn{conn1}.roi1 == pconn{conn2}.roi1 || pconn{conn1}.roi2 == pconn{conn2}.roi2)
                val = 1;
            end
            
            if (pconn{conn1}.roi1 == pconn{conn2}.roi2 || pconn{conn1}.roi2 == pconn{conn2}.roi1)
                val = 1;
            end
            
            % store output
            out{ii} = val;
            
        otherwise
            
            error('Invalid metric between edges requested.');
            
    end
    
end
time = toc;

disp(['Computed link network metrics in ' num2str(time/60) ' minutes.' ]);

disp('Creating link network graph...');

% preallocate matrix
omat = nan(lnodes, lnodes, nmiss);

% create outputs
for ii = 1:size(indx, 1)
    
    % grab the precomputed indices
    conn1 = indx(ii, 1);
    conn2 = indx(ii, 2);
    
    for jj = 1:nmiss
        
        % create matrix
        omat(conn1, conn2, jj) = out{ii}(jj);
        omat(conn2, conn1, jj) = out{ii}(jj);
    
    end
    
end

% fix nan/inf/neg values to zero
omat = squeeze(omat);
omat(isinf(omat)) = 0;
omat(isnan(omat)) = 0;
%omat(omat < 0) = 0;

end

