function [ netw, out ] = fnCreateTreeEdges(parc, fg, maxd, names, minNum, varargin)
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
% Brent McPherson (c), 2021 - Indiana University
%

%% parse arguments

% grab name / value pairs for optional streamline measures
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

% check if fibers are shaped correctly
if size(fg.fibers{1}, 1) ~= 3
    error('Expecting fibers with size(3, N)'); 
end

% set a default maximum distance in mm an endpoint can be from a node.
if(~exist('maxd', 'var') || isempty(maxd))
    maxd = 4; % mrtrix3 uses 4
end

% by default don't enforce a minimum number of streamlines for an edge to get stored
if(~exist('minNum', 'var') || isempty(minNum))
    minNum = 0;
end

%% extract fiber endpoints and parcellation points in acpc space

disp('Converting parcellation labels to AC-PC space...');

% catch xform matrices from parcellation
parc_acpc2img = niftiGet(parc, 'qto_ijk');
parc_img2acpc = niftiGet(parc, 'qto_xyz');

% grab voxel resolution
dvoxmm = abs(parc_img2acpc([1, 6, 11]));

% pull the dimensions of the image
dsize = size(parc.data);

% pull all non-zero labels
labels = unique(parc.data);
labels = labels(labels > 0);
nlabels = size(labels, 1);

% check against what labels should exist in case a label is missing?
% - that's a pain and should be caught before this point, frankly...

% check if names exist / make empty of they don't
if(~exist('names', 'var') || isempty(names))
    disp('No ROI names passed. Will fill in empty names...');
    names = cell(size(labels));
end

% sanity check if the number of names passed equals the number of labels
if(length(names) ~= nlabels)
    warning('The number of provided names do not match the number of labels in parc. Provided names will not be stored.');
    names = cell(size(labels));
end

% pull all the non-zero label voxels from parcellation
pidx = parc.data(:) > 0;

% grab 'gm' volume based on parcellation dimensions
gmvol = sum(pidx) * prod(dvoxmm);

% pull all the non-zero parcellation values, same order as the non-zero indices
pval = parc.data(pidx);

% turn the indices into the cooresponding image coordinates
[ i, j, k ] = ind2sub(dsize, find(pidx));

% convert i/j/k coordiates of parcellation data to AC-PC space
pdat = mrAnatXformCoords(parc_img2acpc, [ i j k ]);

clear pidx i j k

% create the kdtree for assignment of endpoints to labels
ns = createns(pdat, 'nsmethod', 'kdtree');

disp('Computing length of every streamline in fg...')
fblen = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), fg.fibers, 'UniformOutput', true);

disp('Assigning streamlines to parcellation labels...');

% pull all endpoints from fg in an elegant enough way
ep = cell2mat(cellfun(@(x) [ x(:,1)' x(:,end)' ], fg.fibers, 'UniformOutput', false));
ep = [ ep(:,1:3); ep(:,4:6) ];

% get a cellarray of sorted label indices w/in mdist of endpoints
idx = rangesearch(ns, ep, maxd, 'SortIndices', true);

% fill empty indices with a single 0 - makes it easier to iterate over
idx(cellfun(@isempty, idx)) = {0};

% pull the endpoint values and increment by 1 
epi = cellfun(@(x) x(1), idx) + 1; % zero (empty indices) become 1

% add zero value for unassigned endpoints (index 1 becomes an empty code)
epvl = [ 0; pval ]; 

% grab all endpoint labels, 0 is empty
epv = epvl(epi);

% split back into first/second endpoints in fg order
epv1 = epv(1:nfib);
epv2 = epv(nfib+1:end);

clear fib idx ep ep1 ep2 epi epvl epv
% ep/epv useful for inspecting the assignment plot

% % display assigned / unassigned streamlines next to parc in ac-pc space
% epa = epv ~= 0; epm = epv == 0;
% figure; hold on;
% plot3(pdat(:,1), pdat(:,2), pdat(:,3), '.', 'Color', [ 0 0 0.75 ], 'MarkerSize', 3);
% plot3(ep(epa,1), ep(epa,2), ep(epa,3), '+', 'Color', [ 0 0.75 0 ], 'MarkerSize', 2);
% plot3(ep(epm,1), ep(epm,2), ep(epm,3), '+', 'Color', [ 0.75 0 0 ], 'MarkerSize', 2);
% axis image; set(gca, 'Zlim', [ 10 15 ]);

%% store in the network object for reference

netw.parc.xform.acpc2img = parc_acpc2img;
netw.parc.xform.img2acpc = parc_img2acpc;
netw.parc.dsize = dsize;
netw.parc.voxmm = dvoxmm;
netw.parc.labels = labels;
netw.volume.parc = gmvol;

%% catch data about nodes

% preallocate outputs
nodes = cell(length(labels), 1);
tfib = zeros(length(labels), 1);
dfib = zeros(length(labels), 1);
dcnt = zeros(length(labels), 1);

% for every node, assign values volume and duplicate streamline info
for lab = 1:nlabels
    
    % catch info about the labels
    nodes{lab}.name = names(lab);
    nodes{lab}.label = labels(lab);
    nodes{lab}.size = sum(pval == labels(lab));
    nodes{lab}.volume = nodes{lab}.size * prod(dvoxmm);
    nodes{lab}.prop = nodes{lab}.volume / gmvol;
    
    % grab logical of index assignments
    epl1 = epv1 == labels(lab);
    epl2 = epv2 == labels(lab);
    
    % find when both endpoints of a streamline are in the node
    beps = find(epl1 & epl2);
    nodes{lab}.botheps.indices = beps;
    nodes{lab}.botheps.length = fblen(beps);
    
    % catch all values for streamlines with both eps
    for jj = 1:size(nam, 2) 
        nodes{lab}.botheps.(nam{jj}) = val{jj}(beps);
    end
    dfib(lab) = size(beps, 1);
    
    % keep track of labels with both end points of streamlines 
    if ~isempty(beps)
        dcnt(lab) = labels(lab);
    end
        
    % create ROI center in image and acpc space
    center = mean(pdat(pval == labels(lab),:));
    nodes{lab}.center.img = round(mrAnatXformCoords(parc_acpc2img, center)) + 1;
    nodes{lab}.center.acpc = center;
    
    % if there are not any endpoints in this label, throw a warning
    if ~any(epl1 | epl2)
        warning(['ROI label ' num2str(labels(lab)) ' has no streamline terminations.']);
    end
    
    % total fibers assigned to an endpoint
    tfib(lab) = sum(epl1) + sum(epl2) - dfib(lab);
    
end

disp([ 'Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*nfib) ' terminations.' ]);
disp([ num2str(length(dcnt)) ' ROIs had both terminations of ' num2str(sum(dfib)) ' total streamlines.']);

clear ii jj epl1 epl2 beps center

% assign nodes array to output
netw.nodes = nodes;

%% build network edges

% build every unique combination of labels - store preallocated indices of upper diagonal
pairs = nchoosek(1:length(labels), 2);
netw.parc.pairs = pairs;
nedges = size(pairs, 1);

% preallocate paired connection object
edges = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);
ncon = zeros(length(pairs), 1);

disp('Building network edges...');

% for every pair of nodes, estimate the connection
for edge = 1:nedges
    
    % create shortcut names
    roi1 = labels(pairs(edge, 1));
    roi2 = labels(pairs(edge, 2));
        
    % create logical indices for endpoints
    epe1 = epv1 == roi1 & epv2 == roi2;
    epe2 = epv1 == roi2 & epv2 == roi1;
    
    % assign intersections of terminating streamlines
    edges{edge}.fibers.indices = find(epe1 | epe2);
    edges{edge}.fibers.lengths = fblen(edges{edge}.fibers.indices);
    
    % for every variable input, assign it to the connection
    for jj = 1:size(nam, 2)
        edges{edge}.fibers.(nam{jj}) = val{jj}(edges{edge}.fibers.indices);
    end
    
    % enforce a minimum number of streamlines in an edge from assignment
    if size(edges{edge}.fibers.indices, 1) <= minNum
        edges{edge}.fibers = structfun(@(x) [], edges{edge}.fibers, 'UniformOutput', false);
    end
        
    % keep running total of streamlines assigned to a connection
    tcon(edge) = size(edges{edge}.fibers.indices, 1);
    
    % if the connection is empty or not
    if size(edges{edge}.fibers.indices, 1) > 0
        ncon(edge) = 1;
    end
    
end

clear ii jj roi1 roi2 epe1 epe2

% store edges in network
netw.edges = edges;

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
out.conn.possible = size(edges, 1);

% within node connection summary
out.dupf.count = sum(dfib);
out.dupf.roi = dcnt;

end
