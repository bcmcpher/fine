function [ netw, out ] = fnCreateEdges(parc, fg, names, minNum, varargin)
%feCreatePairedConnections creates pconn object of every possible unique pair of 
% of labels in a parcellation w/ streamlines from fibers object.
%
% INPUTS:
%     parc      - loaded, labeled nifti object; the nodes of the network
%     fibers    - cell array of the acpc transformed fibers field of a fiber group
%     fibLength - vector of the length of each streamline
%     weights   - vector of the weights from each streamline 
%
% OUTPUTS:
%     pconn - paired connections cell array, contains data for each edge
%     rois  - data computed for each node of the network
%     out   - catch simple counts of assigned streamlines, connections, etc.
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
% Brent McPherson (c), 2017 - Indiana University
%

%% parse arguments

% grab name / value pairs
nam = varargin(1:2:end);
val = varargin(2:2:end);

% make sure all names are strings
if ~all(cellfun(@(x) ischar(x), nam))
    error('Optional data should be passed in (''name'', value) order.');
end

% make sure all data arrays are numeric
if ~all(cellfun(@(x) isnumeric(x), val))
    error('Optional data arrays must be numeric values.');
end

% make sure all data arrays are the right size
nfib = size(fg.fibers, 1);
if ~all(cellfun(@(x) size(x, 1) == nfib, val))
    error('Optional data arrays must all have x dimension %d.', nfib);
end

% by default don't enforce a minimum number of streamlines for an edge to get stored
if(~exist('minNum', 'var') || isempty(minNum))
    minNum = 0;
end

%% extract fibers to acpc space and identify endpoint coordinates

disp('Computing length of every streamline in fg...')
fibLength = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), fg.fibers, 'UniformOutput', true);

disp('Assigning streamline endpoints to ROI labels...')

% catch xform matrices for parc - add to master out
parc_acpc2img = niftiGet(parc, 'qto_ijk');
parc_img2acpc = niftiGet(parc, 'qto_xyz');

% grab voxel resolution - add to master out
dvoxmm = abs(parc_img2acpc([1, 6, 11]));

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

% check if names exist / make empty of they don't
if(~exist('names', 'var') || isempty(names))
    disp('No ROI names passed. Will fill in empty names...');
    names = cell(size(labels, 1));
end

disp(['Matching streamlines to ' num2str(length(labels)) ' nodes...']);

% preallocate outputs
rois = cell(length(labels), 1);
indx = cell(length(labels), 1);
tfib = zeros(length(labels), 1);
dfib = zeros(length(labels), 1);
dcnt = zeros(length(labels), 1);

% grab size and data of labels
parc_data = parc.data;
gmvol = sum(parc.data(:) > 0) * prod(dvoxmm);

% store in the network object for reference
netw.parc.xform.acpc2img = parc_acpc2img;
netw.parc.xform.img2acpc = parc_img2acpc;
netw.parc.dsize = size(parc.data);
netw.parc.voxmm = dvoxmm;
netw.parc.labels = labels;
netw.volume.parc = gmvol;

% for every label, assign endpoints
for ii = 1:length(labels)
    
    % catch label info
    rois{ii}.name = names(ii);
    rois{ii}.label = labels(ii);
    
    % pull indices for a label in image space
    [ x, y, z ] = ind2sub(netw.parc.dsize, find(parc_data == labels(ii)));
    imgCoords   = [ x, y, z ];

    % catch size of ROI
    rois{ii}.size = size(unique(imgCoords, 'rows'), 1);
    rois{ii}.volume = rois{ii}.size * prod(dvoxmm);
    rois{ii}.prop = rois{ii}.volume / gmvol;
    
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
    bep = intersect(roi_iep1, roi_iep2);
    rois{ii}.botheps.indices = bep;
    rois{ii}.botheps.length = fibLength(bep);
    for jj = 1:size(nam, 2) % catch all values for streamlines with both eps
        rois{ii}.botheps.(nam{jj}) = val{jj}(bep);
    end
    dfib(ii) = size(bep, 1);
    
    % keep track of labels with both end points of streamlines 
    if ~isempty(bep)
        dcnt(ii) = labels(ii);
    end
        
    % create ROI centroid
    rois{ii}.centroid.img = round(mean(imgCoords, 1) + 1);
    rois{ii}.centroid.acpc = mrAnatXformCoords(parc_img2acpc, rois{ii}.centroid.img);
    
    % throw a warning if no terminations are in a label
    if isempty(indx{ii}.indices)
        warning(['ROI label ' num2str(labels(ii)) ' has no streamline terminations.']);
    end
    
    % total fibers assigned to an endpoint
    tfib(ii) = length(indx{ii}.indices);
    
end

disp([ 'Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*nfib) ' terminations.' ]);
disp([ num2str(length(dcnt)) ' ROIs had both terminations of ' num2str(sum(dfib)) ' total streamlines.']);

clear ii x y z imgCoords roi_ep1 roi_ep2 fibers time

netw.nodes = rois;

%% build paired connections object

% build every unique combination of labels
pairs = nchoosek(1:length(labels), 2);

% preallocate paired connection object
pconn = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);
ncon = zeros(length(pairs), 1);

disp('Building paired connections...');

% keep preallocated indices of upper diagonal
netw.parc.pairs = pairs;

% for every pair of nodes, estimate the connection
for ii = 1:length(pairs)
    
    % create shortcut names
    roi1 = indx{pairs(ii, 1)}.indices;
    roi2 = indx{pairs(ii, 2)}.indices;
        
    % assign intersections of terminating streamlines
    pconn{ii}.fibers.indices = intersect(roi1, roi2);
    pconn{ii}.fibers.lengths = fibLength(pconn{ii}.fibers.indices);
    
    % for every variable input, assign it to the connection
    for jj = 1:size(nam, 2)
        pconn{ii}.fibers.(nam{jj}) = val{jj}(pconn{ii}.fibers.indices);
    end
    
    % enforce a minimum number of streamlines in an edge from assignment
    if size(pconn{ii}.fibers.indices, 1) < minNum
        pconn{ii}.fibers = structfun(@(x) [], pconn{ii}.fibers, 'UniformOutput', false);
    end
        
    % keep running total of streamlines assigned to a connection
    tcon(ii) = size(pconn{ii}.fibers.indices, 1);
    
    % if the connection is empty or not
    if size(pconn{ii}.fibers.indices, 1) > 0
        ncon(ii) = 1;
    end
    
end

netw.edges = pconn;

disp(['Built network object with ' num2str(sum(ncon)) ' edges containing ' num2str(sum(tcon)) ' streamlines.']);

%% simple summary of counts

% streamline summaries
out.streamlines.edges = sum(tcon);
out.streamlines.total = nfib;

% end point summaries
out.ep.total = sum(tfib);
out.ep.possible = out.ep.total / (nfib * 2);

% connection summaries
out.conn.total = sum(ncon);
out.conn.possible = size(pconn, 1);

% within node connection summary
out.dupf.count = sum(dfib);
out.dupf.roi = dcnt;

end
