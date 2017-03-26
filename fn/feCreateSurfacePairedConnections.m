%function [ pconn, rois ] = feCreateSurfacePairedConnections(fe, aparc)
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
%[ lcoord, lface ] = read_surf('test/lh.white');
%[ rcoord, rface ] = read_surf('test/rh.white');

% load labels - vertex index, region label, labeled color structure
%[ lvert, ltab, lcol ] = read_annotation('test/lh.aparc.DKTatlas40.annot');
%[ rvert, rtab, rcol ] = read_annotation('test/rh.aparc.DKTatlas40.annot');

% NOW GO FROM HCP SURFS / GLASSER SURF LABELS
% THEY CATCH THE LABELS, OH MY GOD!!!!

% from gifti, import geometry and labels
lhshape = gifti('test/lh.white.surf.gii');
rhshape = gifti('test/rh.white.surf.gii');

lhlabel = gifti('test/lh.aparc.DKTatlas40.label.gii');
rhlabel = gifti('test/rh.aparc.DKTatlas40.label.gii');

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

% merge endpoints into one large vector
ep = [ ep1; ep2 ];

clear ii ep1 ep2

%% assign fiber endpoints to labels

% find number of streamlines
nep = size(fg.fibers, 1);

% combine left and right coordinates
surf_coord = [ lhshape.vertices; rhshape.vertices ];

% combine and make sure left/right labels are separated
%surf_label = [ lhlabel.cdata + 999; rhlabel.cdata + 1999 ]; % HCP data didn't like this
surf_label = [ lhlabel.cdata; rhlabel.cdata ];

% define search radius around streamline termination
radius = 4;

% define bounding box around streamline endpoints that will be searched
bbox = 6;

% preallocate output for parallelized run
out = zeros(nep * 2, 1);

% for every streamline endpoint
tic;
parfor ii = 1:size(ep, 1)
    
    % catch streamline termination
    pt = ep(ii, :);
    
    % create bounding box of + / - bbox mm for faster search
    pt_upper = pt + bbox;
    pt_lower = pt - bbox;
    
    % apply bounding box to search space - creates logical index of all
    % vertices in the surfaces that can be sampled for labels.
    surf_logic = surf_coord(:, 1) < pt_upper(1) & surf_coord(:, 1) > pt_lower(1) & ...
                 surf_coord(:, 2) < pt_upper(2) & surf_coord(:, 2) > pt_lower(2) & ...
                 surf_coord(:, 3) < pt_upper(3) & surf_coord(:, 3) > pt_lower(3);
    
    % if there's nothing in the bounding box, move on
    if (sum(surf_logic) == 0)
        continue
    end
             
    % create bound box of coordinates and matching labels
    tmp_space = surf_coord(surf_logic, :);
    tmp_label = surf_label(surf_logic);
    
    % find the tmp vertices within the space
    tmp_srchs = rangesearch(tmp_space, ep(ii, :), radius);
    
    % if there's nothing within the search radius, move on
    if isempty(tmp_srchs{:})
        out(ii) = nan;
        continue
    end
    
    % find the mode of the found vertices labels as the streamline termination
    lab = mode(tmp_label(tmp_srchs{:}));
    
    % assign label to output
    out(ii) = lab;
    
end
time = toc;

display(['Successfully assigned labels from the surface for ' num2str(2*size(fe.life.M.Phi, 3)) ' endpoints in ' num2str(round(time)/60) ' minutes.']);

% recreate an indexed structure so streamline indices can be assigned
epout = reshape(out, [nep 2]);

% create the unique list of labels on the surface
labels = unique(surf_label);

% create rois cell array for each label in surface vertices
rois = cell(length(labels), 1);
tfib = zeros(length(labels), 1);

display('Sorting assigned streamlines by surface labels...');

tic;
for ii = 1:size(labels, 1)
    
    % catch label info
    rois{ii}.label = labels(ii);
    
    % what do I do about ROI size?
    % for now it is the number of vertices - breaks density estimate
    rois{ii}.size = sum(surf_label == labels(ii));
    
    % assign labels to endpoint indices
    x1 = find(epout(:, 1) == labels(ii));
    x2 = find(epout(:, 2) == labels(ii));
    
    % capture properties I'm interested in
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

display(['Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*size(fe.life.M.Phi, 3)) ' endpoints in ' num2str(round(time)/60) ' minutes.']);

%% CHECK THAT THIS LOOP MATCHES THE VOLUME VERSION

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
    pconn{ii}.all.indices = intersect(roi1.indices, roi2.indices);
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

%end
