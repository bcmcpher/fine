function [ pconn, rois ] = feCreateSurfacePairedConnections(fe, aparc)
%feCreatePairedConnections creates pconn object of every possible unique pair of 
% of labels on a labeled cortical surface. 
%   
%   NOT WORKING - RUN AS SCRIPT
%
%   TODO:
%   - actually have to solve this... 
%

%% load data

% if string is passed, assume it's a path and load it
if isstring(fe) 
    display('Loading fe data...');
    load(fe);
end

% load surfaces - x/y/z vertex in acpc, face data
[ lcoord, lface ] = read_surf('test/lh.white');
[ rcoord, rface ] = read_surf('test/rh.white');

% load labels - vertex index, region label, labeled color structure
[ lvert, ltab, lcol ] = read_annotation('test/lh.aparc.DKTatlas40.annot');
[ rvert, rtab, rcol ] = read_annotation('test/rh.aparc.DKTatlas40.annot');

%% extract fibers to acpc space and identify endpoint coordinates

display('Coverting streamlines and ROIs to ACPC space...')

% convert fibers to acpc space
fg = feGet(fe, 'fg acpc');

% get fiber lengths
fibLength = fefgGet(fg, 'length');

% initialize endpoint outputs
ep1 = zeros(length(fg.fibers), 3);
ep2 = zeros(length(fg.fibers), 3);

% for every fiber, pull the end points
for ii = 1:length(fg.fibers)
    ep1(ii,:) = fg.fibers{ii}(:,1)';
    ep2(ii,:) = fg.fibers{ii}(:,end)';
end

clear ii

%% assign fiber endpoints to labels

% find number of end points
nep = length(fg.fibers);

% combine left and right coordinates / labels
% make sure left/right labels are separated
surf_coord = [ lcoord; rcoord ];
surf_label = [ ltab; rtab + 1 ];

% define search radius around streamline termination
radius = 3;

% combine unique label indices
%labels = [ unique(ltab(ltab > 0)); unique(rtab(rtab > 0)) + 1 ];

% preallocate output for parallelized run
out = cell(nep, 1);

% by endpoint
parfor ii = 1:nep
    
    % catch first terminations
    pt1 = ep1(ii, :);
    pt2 = ep2(ii, :);
    
    % create bounding box of + / - 5mm for faster search
    pt1_upper = pt1 + 5;
    pt1_lower = pt1 - 5;
    
    pt2_upper = pt2 + 5;
    pt2_lower = pt2 - 5;
    
    % apply bounding box to search space
    % DO IN 1 STATEMENT FOR SINGLE LOGICAL INDEX - NEED TO NOT BREAK LABELS
    surf_space1 = surf_coord(surf_coord(:, 1) < pt1_upper(1) & surf_coord(:, 1) > pt1_lower(1), :);
    surf_space1 = surf_space1(surf_space1(:, 2) < pt1_upper(2) & surf_space1(:, 2) > pt1_lower(2), :);
    surf_space1 = surf_space1(surf_space1(:, 3) < pt1_upper(3) & surf_space1(:, 3) > pt1_lower(3), :);
    
    surf_space2 = surf_coord(surf_coord(:, 1) < pt2_upper(1) & surf_coord(:, 1) > pt2_lower(1), :);
    surf_space2 = surf_space2(surf_space2(:, 2) < pt2_upper(2) & surf_space2(:, 2) > pt2_lower(2), :);
    surf_space2 = surf_space2(surf_space2(:, 3) < pt2_upper(3) & surf_space2(:, 3) > pt2_lower(3), :);
    
    % THIS BREAKS THE INDICES OF THE LABELS - HAVE TO PERFORM EQUIVALENT
    % SUBSET OR RE-INDEX INTO THE ENTIRE LIST
    
%     % find vertices within radius of cortex node 
%     out{ii}.rsrch1 = rangesearch(surf_coord, ep1(ii, :), radius);
%     out{ii}.rsrch2 = rangesearch(surf_coord, ep2(ii, :), radius);

    % find vertices within bounding box of space
    out{ii}.rsrch1 = rangesearch(surf_space1, ep1(ii, :), radius);
    out{ii}.rsrch2 = rangesearch(surf_space2, ep2(ii, :), radius);
    
    
    % for end points with vertices, find the mode of vertex labels
    if ~isempty(out{ii}.rsrch1{:})
        lab1 = mode(surf_label(out{ii}.rsrch1{:}));
    else
        lab1 = nan;
    end
    
    if ~isempty(out{ii}.rsrch2{:})
        lab2 = mode(surf_label(out{ii}.rsrch2{:}));
    else
        lab2 = nan;
    end
    
    % catch the label in a matrix to create the roi summary I want
    out{ii}.label1 = lab1;
    out{ii}.label2 = lab2;
    
end
    
epout = zeros(nep, 2);
% create matrix of all labels for fast creation
for ii = 1:length(out)
    epout(ii, 1) = out{ii}.label1; 
    epout(ii, 2) = out{ii}.label2;
end

% create rois cell array for each label in surface vertices
rois = cell(length(labels), 1);
tfib = zeros(length(labels), 1);

tic;
for ii = 1:length(labels)
    
    % catch label info
    rois{ii}.label = labels(ii);
    
    % assign labels to endpoint indices
    x1 = find(epout(:, 1) == labels(ii));
    x2 = find(epout(:, 2) == labels(ii));
    
    % campture properties I'm interested in
    rois{ii}.indices = unique([ x1; x2 ]);
    rois{ii}.lengths = fibLength(rois{ii}.indices);
    rois{ii}.weights = fe.life.fit.weights(rois{ii}.indices);
    
    if isempty(rois{ii}.indices)
        warning(['ROI label ' num2str(labels(ii)) ' is not assigned any streamline terminations.']);
    end
    
    % add roi surface centroids
    % RECONCILE THAT I NEED TO MANIPULATE L/R LABELS
    % MAKE SURE ACPC2IMG IS IMPLEMENTED CORRECTLY - PROBABLY NEED IMAGE HEADER
    %rois{ii}.centroid.acpc = mean(surf_coords(surf_label == label(ii)));
    %rois{ii}.centroid.img = fe.life.xform.acpc2img * rois{ii}.centroid.acpc;
    
    % counter of how many streamline terminations are assigned
    tfib(ii) = size(rois{ii}.indices, 1);
    
end
time = toc;

display(['Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*size(fe.life.M.Phi, 3)) ' endpoints in ' num2str(round(time)) ' seconds.']);
% no rh rois and 4 lh rois are not assigned terminations despite assigning
% over 420000 streamlines. Figure out why after maintenance.

% % for every unique label, find endpoint vertices within minimum distance of labeled surface vertex
% for ii = 1:length(labels)
%     
%     % pull the vertices for the label
%     labVerts = lcoord(ltab == labels(ii), :);
%     labIndex = find(ltab == labels(ii));
%     
%     % preallocate vertice catcher
%     circaVert1 = cell(nep, 1);
%     circaVert2 = cell(nep, 1);
%     
%     % check for each end point if it's within distance
%     parfor jj = 1:nep
%         circaVert1{jj} = rangesearch(labVerts, ep1(jj, :), 3);
%         circaVert2{jj} = rangesearch(labVerts, ep2(jj, :), 3);
%     end
%     
%     % find indices that are not empty
%     for jj = 1:nep
%         circaIndx1(jj) = ~isempty(circaVert1{jj, 1}{:});
%         circaIndx2(jj) = ~isempty(circaVert2{jj, 1}{:}); 
%     end
%     
%     % how many are found - reasonable...
%     %sum(circaIndx1) + sum(circaIndx2)
%     
%     % catch indices
%     out{ii}.indices = [ find(circaIndx1) find(circaIndx2) ];
%     
%     % THIS DOES NOT STOP DOUBLE COUNTING
%     
% end

% % ORIGINAL END POINT ASSIGNMENT
% display(['Matching streamlines to ' num2str(length(labels)) ' nodes...']);
% 
% % preallocate outputs
% rois = cell(length(labels), 1);
% tfib = zeros(length(labels), 1);
% 
% % for every label, assign endpoints
% tic;
% for ii = 1:length(labels)
%     
%     % catch label info
%     rois{ii}.label = labels(ii);
%     
%     % pull indices for a label in image space
%     [ x, y, z ] = ind2sub(size(aparc.data), find(aparc.data == labels(ii)));
%     imgCoords = [ x, y, z ];
%     
%     % convert label indices to ACPC coordinates
%     acpcCoords = mrAnatXformCoords(aparc_img2acpc, imgCoords); % rely on file header
%     acpcCoords = round(acpcCoords) + 1;
%    
%     % catch size 
%     rois{ii}.size = size(unique(acpcCoords, 'rows'), 1);
%     
%     % find streamline endpoints in ROI acpc coordinates
%     roi_ep1 = ismember(ep1, acpcCoords, 'rows');
%     roi_ep2 = ismember(ep2, acpcCoords, 'rows');
%     
%     % for fibers that end in rois, catch indices / lengths / weights 
%     fibers = [ find(roi_ep1); find(roi_ep2) ];
%     rois{ii}.end.fibers = fibers;
%     rois{ii}.end.length = fibLength(rois{ii}.end.fibers);
%     rois{ii}.end.weight = fe.life.fit.weights(rois{ii}.end.fibers);
%         
%     % create ROI centroid
%     rois{ii}.centroid.acpc = round(mean(acpcCoords) + 1);
%     rois{ii}.centroid.img = round(mean(imgCoords)) + 1;
%     
%     % create endpoint density ROI object
%     
%     if isempty(rois{ii}.end.fibers)
%         warning(['ROI label ' num2str(labels(ii)) ' has no streamline terminations.']);
%     end
%     
%     % total fibers assigned to an endpoint
%     tfib(ii) = length(rois{ii}.end.fibers);
%     
% end
% time = toc; 
% 
% display(['Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*size(fe.life.M.Phi, 3)) ' endpoints in ' num2str(round(time)) ' seconds.']);
% 
% clear ii x y z imgCoords acpcCoords roi_ep1 roi_ep2 fibers tfib time

%% build  paired connections object

% build every unique combination of labels
pairs = nchoosek(1:length(labels), 2);

% preallocate paired connection object
pconn = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);

display('Building paired connections...');

% for every pair of nodes, estimate the connection
for ii = 1:length(pairs)
    
    % create shortcut names
    roi1 = rois{pairs(ii, 1)};
    roi2 = rois{pairs(ii, 2)};
    
    % catch shortcut names
    pconn{ii}.roi1 = roi1.label;
    pconn{ii}.roi2 = roi2.label;
    
    % catch size of ROIs
    pconn{ii}.roi1sz = roi1.size;
    pconn{ii}.roi2sz = roi2.size;
    
    % assign intersections of terminating streamlines
    pconn{ii}.all.indices = intersect(roi1.end.fibers, roi2.end.fibers);
    pconn{ii}.all.lengths = fibLength(pconn{ii}.all.indices);
    pconn{ii}.all.weights = fe.life.fit.weights(pconn{ii}.all.indices);
    
    % find all weighted fibers
    nzw = pconn{ii}.all.weights > 0;
    
    pconn{ii}.nzw.indices = pconn{ii}.all.indices(nzw);
    pconn{ii}.nzw.lengths = pconn{ii}.all.lengths(nzw);
    pconn{ii}.nzw.weights = pconn{ii}.all.weights(nzw);
    
    % calculate combined size of ROI
    psz = pconn{ii}.roi1sz + pconn{ii}.roi2sz;
    
    % values used to calculate different edge weights
    cnt = size(pconn{ii}.all.indices, 1);
    len = mean(pconn{ii}.all.lengths);
    dln = sum(1 / pconn{ii}.all.lengths);
    
    % values used to calculate different edge weights for non-zero weighted streamlines
    nzcnt = size(pconn{ii}.all.indices(nzw), 1);
    nzlen = mean(pconn{ii}.all.lengths(nzw));
    if isempty(pconn{ii}.nzw.lengths) % if there are no nz lengths
        nzdln = 0;
    else
        nzdln = sum(1 / pconn{ii}.nzw.lengths);
    end
    
    % create all streamline counts
    pconn{ii}.all.matrix.count = cnt;
    pconn{ii}.all.matrix.density = (2 * cnt) / psz;
    pconn{ii}.all.matrix.length = len;
    pconn{ii}.all.matrix.denlen = (2 / psz) * dln;
        
    % create non-zero weighted streamline counts
    pconn{ii}.nzw.matrix.count = nzcnt;
    pconn{ii}.nzw.matrix.density = (2 * nzcnt) / psz;
    pconn{ii}.nzw.matrix.length = nzlen;
    pconn{ii}.nzw.matrix.denlen = (2 / psz) * nzdln;
    
    % keep running total of total streamlines assigned a connection
    tcon(ii) = size(pconn{ii}.all.indices, 1);
    
end

display(['Built paired connections object with ' num2str(sum(tcon)) ' streamlines.']);

end
