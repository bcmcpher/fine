function [ pfg ] = feCleanConnectomeEdges(parc, fg, minLength, maxDist, maxLengthStd, numNodes, maxIter, minNum)
%
% STILL IN DEVELOPMENT - THE COUNT IT PRINTS DOES NOT MATCH WHAT IT DOES
%
%feCleanPairedConnections cleans the outlier fibers from connections inside 
% a paired connection (pconn) cell array. 
%
% INPUTS:
%    fg        - fiber group in acpc space
%    pconn     - paired connection object created with fg
%    label     - string indicating the fiber groups to clean
%                either:
%                       'all' for all assigned streamlines or
%                       'nzw' for non-zero weighted fibers returned by LiFE
%    minLength    - the minimum length in mm of a connection to be considered
%                   valid (default = 10)
%    maxDist      - the # of standard deviations a streamline node can take from the
%                   centroid of the tract (default = 4)
%    maxLengthStd - the # of standard deviations the length of a
%                   streamline can deviate from the tracts mean (default = 4)
%    numNodes     - the number of nodes fibers are resampled to for
%                   similarity comparisons (default = 100)
%    maxIter      - the number of iterations outliers will be tested and
%                   sampled from (default = 10)
%    minNum       - the minimum number of streamlines that must exist for
%                   a connection to be kept (default = 3).
%                   NOTE: tract profiles require a minimum of 3 streamlines.
%             
% OUTPUTS:
%     pconn - is the paired connections object with the cleaned streamline
%             indices added in a new label appended with '_clean'.
%
%     cln   - debugging output; the cell array that is added internally to pconn
%
% EXAMPLE:
%
% % load data
% parc = niftiRead('labels.nii.gz');
% fg        = feGet(fe, 'fibers acpc');
% fibers    = fg.fibers;
% fibLength = fefgGet(fg, 'length');
% weights   = feGet(fe, 'fiberweights');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % clean streamlines too far from center with all streamlines that are assigned
% pconn = feCleanPairedConnections(fg, pconn, 'all');
%
% Brent McPherson (c), 2017 - Indiana University
%

% fill in defaults for any missing parameters

if(~exist('minLength', 'var') || isempty(minLength))
    minLength = 10;
end

if(~exist('maxDist', 'var') || isempty(maxDist))
    maxDist = 4;
end

if(~exist('maxLengthStd', 'var') || isempty(maxLengthStd))
    maxLengthStd = 4;
end

if(~exist('numNodes', 'var') || isempty(numNodes))
    numNodes = 100;
end

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('minNum', 'var') || isempty(minNum))
    minNum = 3;
end

%% extract fibers to acpc space and identify endpoint coordinates

% pull total count of streamlines
nfib = size(fg.fibers, 1);

disp('Computing length of every streamline in fg...')
fibLength = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), fg.fibers, 'UniformOutput', true);

disp('Assigning streamline endpoints to ROI labels...')

% catch xform matrices for parc - add to master out
parc_acpc2img = niftiGet(parc, 'qto_ijk');
%parc_img2acpc = niftiGet(parc, 'qto_xyz');

% grab voxel resolution - add to master out
%dvoxmm = abs(parc_img2acpc([1, 6, 11]));

% convert acpc fibers to parcellation space
ifg = dtiXformFiberCoords(fg, parc_acpc2img, 'img');

% re-extract label space fiber coordinates
fibers = ifg.fibers;

% initialize endpoint outputs
ep1 = zeros(length(fibers), 3);
ep2 = zeros(length(fibers), 3);

% for every fiber, pull the end points
if size(fibers{1}, 1) ~= 3
    error('Expecting fibers with size(3, N)'); 
end

% extract every streamline end point
for ii = 1:length(fibers)
    ep1(ii,:) = fibers{ii}(:,1)';
    ep2(ii,:) = fibers{ii}(:,end)';
end

clear ii fibers

% round converted end points to match image indices like LiFE
ep1 = round(ep1) + 1;
ep2 = round(ep2) + 1;

%% assign fiber endpoints to labels

% pull all non-zero labels
labels = unique(parc.data);
labels = labels(labels > 0);

disp(['Matching streamlines to ' num2str(length(labels)) ' nodes...']);

% preallocate outputs
%rois = cell(length(labels), 1);
indx = cell(length(labels), 1);
tfib = zeros(length(labels), 1);
%dfib = zeros(length(labels), 1);
%dcnt = zeros(length(labels), 1);

% grab size and data of labels
parc_data = parc.data;
%gmvol = sum(parc.data(:) > 0) * prod(dvoxmm);

% store in the network object for reference
%netw.parc.xform.acpc2img = parc_acpc2img;
%netw.parc.xform.img2acpc = parc_img2acpc;
dsize = size(parc.data);
%netw.parc.voxmm = dvoxmm;
%netw.volume.gmvol = gmvol;

% for every label, assign endpoints
for ii = 1:length(labels)
    
    % catch label info
    %rois{ii}.name = names(ii);
    %rois{ii}.label = labels(ii);
    
    % pull indices for a label in image space
    [ x, y, z ] = ind2sub(dsize, find(parc_data == labels(ii)));
    imgCoords   = [ x, y, z ];

    % catch size of ROI
    %rois{ii}.size = size(unique(imgCoords, 'rows'), 1);
    %rois{ii}.volume = rois{ii}.size * prod(dvoxmm);
    %rois{ii}.prop = rois{ii}.volume / gmvol;
    
    % find streamline endpoints in image coordinates for label
    roi_ep1 = ismember(ep1, imgCoords, 'rows');
    roi_ep2 = ismember(ep2, imgCoords, 'rows');
    
    % find the indices of the streamlines
    roi_iep1 = find(roi_ep1);
    roi_iep2 = find(roi_ep2);
    
    % combine unique streamline indices
    fibers = unique([ roi_iep1; roi_iep2 ]);
        
    % for fibers that end in rois, catch indices / lengths / weights
    indx{ii}.indices = fibers; % don't store all ROI endoints in ROI
    
    % do both endpoints share the same ROI? count the instances, track indices
    %bep = intersect(roi_iep1, roi_iep2);
    %rois{ii}.botheps.indices = bep;
    %rois{ii}.botheps.length = fibLength(bep);
    %for jj = 1:size(nam, 2) % catch all values for streamlines with both eps
    %    rois{ii}.botheps.(nam{jj}) = val{jj}(bep);
    %end
    %dfib(ii) = size(bep, 1);
    
    % keep track of labels with both end points of streamlines 
    %if ~isempty(bep)
    %    dcnt(ii) = labels(ii);
    %end
    
    % create ROI centroid
    %rois{ii}.centroid.img = round(mean(imgCoords, 1) + 1);
    %rois{ii}.centroid.acpc = mrAnatXformCoords(parc_img2acpc, rois{ii}.centroid.img);
    
    % throw a warning if no terminations are in a label
    if isempty(indx{ii}.indices)
        warning(['ROI label ' num2str(labels(ii)) ' has no streamline terminations.']);
    end
    
    % total fibers assigned to an endpoint
    tfib(ii) = length(indx{ii}.indices);
    
end

disp([ 'Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*nfib) ' terminations.' ]);
%disp([ num2str(length(dcnt)) ' ROIs had both terminations of ' num2str(sum(dfib)) ' total streamlines.']);

clear ii x y z imgCoords roi_ep1 roi_ep2 fibers time

%netw.rois = rois;

%% build paired connections object

% build every unique combination of labels
pairs = nchoosek(1:length(labels), 2);

% preallocate paired connection object
pconn = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);
ncon = zeros(length(pairs), 1);

disp('Building paired connections...');

% keep preallocated indices
%netw.parc.pairs = pairs;

% for every pair of nodes, estimate the connection
for ii = 1:length(pairs)
    
    % create shortcut names
    roi1 = indx{pairs(ii, 1)}.indices;
    roi2 = indx{pairs(ii, 2)}.indices;
    
    % assign intersections of terminating streamlines
    pconn{ii}.fibers.indices = intersect(roi1, roi2);
    pconn{ii}.fibers.lengths = fibLength(pconn{ii}.fibers.indices);
    
    % for every variable input, assign it to the connection
    %for jj = 1:size(nam, 2)
    %    pconn{ii}.fibers.(nam{jj}) = val{jj}(pconn{ii}.fibers.indices);
    %end
    
    % keep running total of total streamlines assigned a connection
    tcon(ii) = size(pconn{ii}.fibers.indices, 1);
    
    % if the connection is empty or not
    if size(pconn{ii}.fibers.indices, 1) > 0
        ncon(ii) = 1;
    end
    
end

%netw.pconn = pconn;

disp(['Built ' num2str(sum(tcon > 0)) ' edges with ' num2str(sum(tcon)) ' streamlines.']);

% KEEP AND PASS THROUGH ANYTHING NOT IN THE NETWORK?
% THIS IS THINGS LIKE CST / CEREBELLUM THAT ARE NOT COMMON IN SIMPLE PARCELLATIONS
% LIKELY NEEDED IF THIS IS A PRE-LiFE STEP

%% simple summary of counts

% % streamline summaries
% out.streamlines.edges = sum(tcon);
% out.streamlines.total = nfib;
% 
% % end point summaries
% out.ep.total = sum(tfib);
% out.ep.possible = out.ep.total / (nfib * 2);
% 
% % connection summaries
% out.conn.total = sum(ncon);
% out.conn.possible = size(pconn, 1);
% 
% % within node connection summary
% out.dupf.count = sum(dfib);
% out.dupf.roi = dcnt;

%% clean the built edges

% add cleaning parameters to output
%netw.clean.minLength = minLength;
%netw.clean.maxDist = maxDist;
%netw.clean.maxLengthStd = maxLengthStd;
%netw.clean.numNodes = numNodes;
%netw.clean.maxIter = maxIter;
%netw.clean.minNum = minNum;

% pull a logical from each edge if any streamlines exist
pc = sum(cellfun(@(x) size(x.fibers.indices, 1) > 0, pconn));
disp(['Started cleaning ' num2str(pc) ' edges...']);

% initialize a count of before / after cleaning edges
bfcln = zeros(size(pconn, 1), 1);
afcln = zeros(size(pconn, 1), 1);

% pull connections from network
%pconn = netw.pconn;

tic;
for ii = 1:size(pconn, 1)
    
    % get the requested edge
    edge = pconn{ii};
    
    % grab the before count
    bfcln(ii) = size(edge.fibers.indices, 1);    
    
    % if the connection isn't empty
    if ~isempty(edge.fibers.indices)
                
        % drop indices less than a particular length
        keep_idx = edge.fibers.lengths > minLength;
                
        % if there are too few streamlines, fill in empty and move on
        if(sum(keep_idx) < minNum)
            edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
            pconn{ii} = edge;
            continue
        else
           % filter other values before fg cleaning
           edge.fibers = structfun(@(x) x(keep_idx), edge.fibers, 'UniformOutput', false);
        end
        
        % create an fg group of the filtered edge
        tfg = fgExtract(fg, edge.fibers.indices, 'keep');
        
        % try and clean the pruned length / count, fail over to empty 
        try
            % for all streamlines, compute outliers
            % cln is a logical of the streamlines to keep
            [ ~, cln ] = mbaComputeFibersOutliers(tfg, maxDist, maxLengthStd, numNodes, 'mean', 0, maxIter);
            
            % filter edge down to cleaned fg
            edge.fibers = structfun(@(x) x(cln), edge.fibers, 'UniformOutput', false);
            
        catch
            warning('Cleaning of edge index %d failed. Assigning empty values in network structure.', ii);
            edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
            pconn{ii} = edge;
            continue
        end
        
        % keep count of cleaned connections
        % (before cleaning) - (total > minLength and after cleaning fxn)
        afcln(ii) = bfcln(ii) - sum(cln);
        
    else 
        
        % fill in empty cells and skip to next connection
        edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
        pconn{ii} = edge;
        continue
        
    end
    
    % just in case it's missed, just fill it back in
    pconn{ii} = edge;
    
end
time = toc;

clear ii edge keep_idx tfg cln

% reassign paired connection
%netw.pconn = pconn;

% compute how many are kept vs. dropped
diffs = bfcln - afcln;

disp(['Cleaned outlier streamlines from ' num2str(sum(afcln ~= bfcln)) ' edges in ' num2str(round(time)/60) ' minutes.']);
disp(['Removed ' num2str(sum(afcln)) ', keeping ' num2str(sum(diffs)) ' streamlines of ' num2str(sum(bfcln)) '.']);

%% prune output fg

% pull a vector of all kept streamlines
ccell = cellfun(@(x) x.fibers.indices, pconn, 'UniformOutput', false);
clean = vertcat(ccell{:});

% copy input fiber group
pfg = fg;

% remove cleaned fibers from fg
pfg.fibers = pfg.fibers(clean);

end
