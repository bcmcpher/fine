function [ netw, out ] = fnCreateEdges(parc, fg, names, maxdist, minNum, varargin)
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
if(~exist('maxdist', 'var') || isempty(maxdist))
    maxdist = 4; % mrtrix3 uses 4
end

% by default don't enforce a minimum number of streamlines for an edge to get stored
if(~exist('minNum', 'var') || isempty(minNum))
    minNum = 0;
end

%% extract fiber endpoints and parcellation points in acpc space

disp('Extracting parcellation information...');

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

%% confirm some correspondance between labels in parc and names

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

%% extract the labels for each entry in parc

% pull all the non-zero label voxels from parcellation
pidx = parc.data(:) > 0;

% grab 'gm' volume based on parcellation dimensions
gmvol = sum(pidx) * prod(dvoxmm);

% pull all the non-zero parcellation values, same order as the non-zero indices
pval = parc.data(pidx);

% turn the indices into the cooresponding image coordinates
[ i, j, k ] = ind2sub(dsize, find(pidx));

% create the voxel coordinate matrix of all pval labels
pijk = [ i, j, k ];
    
% if no distance is allowed, map streamlines to voxels (orig functionality)
if maxdist == 0
    
    disp('Assigning streamline endpoints directly to voxels...');
    
    % convert acpc fibers to parcellation space
    ifg = dtiXformFiberCoords(fg, parc_acpc2img, 'img');
    
    % get endpoints side by side in a [ nfib x 6 ] matrix
    ec = cell2mat(cellfun(@(x) [ x(:,1)' x(:,end)' ], ifg.fibers, 'UniformOutput', false));
    
    % stack and round the transformed endpoints into a single vector
    ep = round([ ec(:,1:3); ec(:,4:6) ]) + 1;
    
    % assign the voxel coordinates for building the kdtree
    pdat = pijk;
    
% otherwise, map parcellation to acpc coords for distance alignment to ep
else
    
    disp(['Assigning streamline endpoints to nearest label within ' num2str(maxdist) 'mm...']);
    
    % get endpoints side by side in a [ nfib x 6 ] matrix
    ec = cell2mat(cellfun(@(x) [ x(:,1)' x(:,end)' ], fg.fibers, 'UniformOutput', false));
    
    % stack the endpoints into a single vector
    ep = [ ec(:,1:3); ec(:,4:6) ];
        
    % convert i/j/k coordiates of parcellation data to AC-PC space
    pdat = mrAnatXformCoords(parc_img2acpc, pijk);
        
end

disp('Assigning streamlines to parcellation labels...');

% create the kdtree for assignment of endpoints to labels
ns = createns(pdat, 'nsmethod', 'kdtree');

% get a cellarray of sorted label indices w/in mdist of endpoints
idx = rangesearch(ns, ep, maxdist, 'SortIndices', true);

% fill empty indices with a single 0 - makes it easier to iterate over
idx(cellfun(@isempty, idx)) = {0};

% pull the endpoint values and increment by 1
epi = cellfun(@(x) x(1), idx) + 1; % zero (empty indices) become 1

% add zero as first index of values for unassigned endpoints 
% index 1 becomes an empty code (0)
epvl = [ 0; pval ];

% grab all endpoint labels, 0 is empty
epv = epvl(epi);

% split back into first/second endpoints in fg order
epv1 = epv(1:nfib);
epv2 = epv(nfib+1:end);

clear pidx i j k fib idx ec ep1 ep2 epi epvl
clear ep epv % ep/epv useful for inspecting the assignment plot below

% % display assigned / unassigned streamlines next to parc in ac-pc space
% epa = epv ~= 0; epm = epv == 0;
% figure; hold on;
% plot3(pdat(:,1), pdat(:,2), pdat(:,3), '.', 'Color', [ 0 0 0.75 ], 'MarkerSize', 3);
% plot3(ep(epa,1), ep(epa,2), ep(epa,3), '+', 'Color', [ 0 0.75 0 ], 'MarkerSize', 2);
% plot3(ep(epm,1), ep(epm,2), ep(epm,3), '+', 'Color', [ 0.75 0 0 ], 'MarkerSize', 2);
% axis image; set(gca, 'Zlim', [ 10 15 ]);

% endpoints when neither label is 0
epu = ~or(epv1 == 0, epv2 == 0);

% endpoints with matching labels
epd = epv1 == epv2; 

% both endoints are assigned to different non-0 labels
ept = epu & ~epd; 

% create sorted ascending order of streamline ep values (upper diagonal)
eps = sort([ epv1, epv2 ], 2);
% the endpoints of only the successfully assigned streamlines

% pull the unique rows of all streamlines that exist
[ uedges, ~, ic ] = unique(eps(ept, :), 'rows');

% get the count for each edge and zero connections below minimum count
h = accumarray(ic, 1); 
h(h <= minNum) = 0; 

% drop unique edge indices w/ a count less than minNum
uedges(h == 0,:) = [];
% uedges is indices into the upper diagonal identifying the edges 
% to be stored  in netw.edges

%assign = [ uedges, h ]; % the uniquely found edges + count

disp('Computing the length of every assigned streamline...')
fblen = zeros(nfib, 1); % use epu b/c ept excluded all self connections, which are desired
fblen(epu) = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), fg.fibers(epu), 'UniformOutput', true);

%% store in the network object for reference

netw.parc.xform.acpc2img = parc_acpc2img;
netw.parc.xform.img2acpc = parc_img2acpc;
netw.parc.dsize = dsize;
netw.parc.voxmm = dvoxmm;
netw.parc.labels = labels;
netw.parc.maxdist = maxdist;
netw.volume.parc = gmvol;

%% catch data about nodes

% preallocate node output and summary counts

% define the empty node array
bept = struct('indices', [], 'length', []);
cntr = struct('img', [], 'acpc', []);

% loop over repeatable input(s) of additional streamline values
for jj = 1:size(nam, 2)
    bept.(nam{jj}) = [];
end

% preallocate paired connection object
nodet = struct('name', [], 'label', [], 'size', [], 'volume', []);
nodet.botheps = bept;
nodet.center = cntr;

% create empty, preallocated cell array of edges
nodes = repmat({nodet}, length(labels), 1);

% iteratively add up what's assigned
tfib = zeros(length(labels), 1);
dfib = zeros(length(labels), 1);
dcnt = zeros(length(labels), 1);

% for every node, assign values volume and duplicate streamline info
for node = 1:nlabels
    
    % catch info about the labels
    nodes{node}.name = names(node);
    nodes{node}.label = labels(node);
    nodes{node}.size = sum(pval == labels(node));
    nodes{node}.volume = nodes{node}.size * prod(dvoxmm);
    
    % find when both endpoints of a streamline are in the node
    beps = find(epd & epv1 == labels(node));
    nodes{node}.botheps.indices = beps;
    nodes{node}.botheps.length = fblen(beps);
    
    % catch all values for streamlines with both eps
    for jj = 1:size(nam, 2) 
        nodes{node}.botheps.(nam{jj}) = val{jj}(beps);
    end
    dfib(node) = size(beps, 1);
    
    % keep track of labels with both end points of streamlines 
    if ~isempty(beps)
        dcnt(node) = labels(node);
    end
        
    % create ROI center in image and acpc space
    center = mean(pijk(pval == labels(node),:));
    nodes{node}.center.img = round(center);
    nodes{node}.center.acpc = mrAnatXformCoords(parc_img2acpc, center);
    
    % if there are not any endpoints in this label, throw a warning
    if ~any(uedges(:) == labels(node))
        warning(['Node label ' num2str(labels(node)) ' has no streamlines assigned.']);
    end
    
    % total fibers assigned to an endpoint
    tfib(node) = sum(eps(:) == labels(node)) - dfib(node);
    
end

disp([ 'Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*nfib) ' streamline endpoints.' ]);
disp([ num2str(length(dcnt)) ' nodes had both endpoints of ' num2str(sum(dfib)) ' total streamlines.']);

clear node jj beps center
clear epl1 epl2 

% assign nodes array to output
netw.nodes = nodes;

%% build network edges

% build every unique combination of labels - store preallocated indices of upper diagonal
nnodes = size(netw.nodes, 1);

% compute indices of output size
mask = triu(ones(nnodes), 1); 

% preallocate unique indices - faster than nchoosek
[ xind, yind ] = find(mask > 0);
pairs = sortrows([ xind, yind ]);

% the total number of possible and non-zero edges to store
nedges = size(pairs, 1);
redges = size(uedges, 1); 

% preallocate edge structure
edget = struct('fibers', struct('indices', [], 'lengths', []));

% loop over repeatable input(s) of additional streamline values
for jj = 1:size(nam, 2)
    edget.fibers.(nam{jj}) = [];
end

% create empty, preallocated cell array of edges
edges = repmat({edget}, redges, 1);

% create zero arrays of counts for output / debugging
tcon = zeros(redges, 1);

disp([ 'Building ' num2str(redges) ' network edges...' ]);

% build label pairs to facilitate skipping empty edges during assignment 
labpr = [ labels(pairs(:,1)) labels(pairs(:,2)) ];

% for every pair of nodes, estimate the connection
for edge = 1:redges
    
    % create shortcut names
    roi1 = uedges(edge, 1);
    roi2 = uedges(edge, 2);
        
    % create logical indices for endpoints
    epe1 = epv1 == roi1 & epv2 == roi2;
    epe2 = epv1 == roi2 & epv2 == roi1;
    
    % the operation to find indices for each nonzero edge
    eidx = find(ept & (epe1 | epe2));

    % assign intersections of terminating streamlines
    edges{edge}.fibers.indices = eidx; %find(epe1 | epe2)
    edges{edge}.fibers.lengths = fblen(eidx);
    
    % for every variable optional input, assign it to the connection
    for jj = 1:size(nam, 2)
        edges{edge}.fibers.(nam{jj}) = val{jj}(eidx);
    end

    % keep running total of streamlines assigned to a connection
    tcon(edge) = size(eidx, 1);
        
end

clear edge jj roi1 roi2 epe1 epe2 eidx

% find the indices of the upper diagonal that are stored in edges
[ ~, idx ] = intersect(labpr, uedges, 'rows');

% store only the indices into the upper diagonal that exist
netw.parc.pairs = pairs(idx,:);

% store edges in network
netw.edges = edges;

disp(['Built network object with ' num2str(redges) ' edges containing ' num2str(sum(tcon)) ' streamlines.']);

%% simple summary of counts

% streamline summaries
out.streamlines.edges = sum(tcon);
out.streamlines.total = nfib;

% end point summaries
out.ep.total = sum(tfib);
out.ep.possible = nfib * 2;

% connection summaries
out.conn.total = redges; %sum(ncon);
out.conn.possible = nedges;

% within node connection summary
out.dupf.count = sum(dfib);
out.dupf.roi = dcnt;

end
