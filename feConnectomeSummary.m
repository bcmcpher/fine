function [ emat, omat, pconn, out ] = feConnectomeSummary(feFile, rois)
%feConnectomeSummary create a network from a fit fe structure and a
%parcellation of cortical regions
%   figure out the rounding
%   streamline what I catch and where
%
% feFile = 'test/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat';
% rois = 'test/aparc+aseg.nii.gz';
%

%% load data

display('Loading data...');

% load the arguments - build checks to use structures if provided
load(feFile);
aparc = niftiRead(rois);

% quick subset of ROIs b/c I don't have the inflated / cleaned ones here
%aparc.data(aparc.data < 1000) = 0;
%aparc.data(aparc.data == 1000 | aparc.data == 2000) = 0;

% catch xform matrices for aparc
aparc_img2acpc = aparc.qto_xyz;
%aparc_acpc2img = aparc.qto_ijk;

%% extract fibers to acpc space and identify endpoint coordinates

% convert fibers to acpc space
%fg = dtiXformFiberCoords(fe.fg, fe.life.xform.img2acpc, 'acpc');
fg = feGet(fe, 'fg acpc');

% initialize endpoint outputs
ep1 = zeros(length(fg.fibers), 3);
ep2 = zeros(length(fg.fibers), 3);

% for every fiber, pull the end points
for ii = 1:length(fg.fibers)
    ep1(ii,:) = fg.fibers{ii}(:,1)';
    ep2(ii,:) = fg.fibers{ii}(:,end)';
end

% round all the endpoints
%ep1 = floor(ep1 + 0.50);
%ep2 = floor(ep2 + 0.50);
ep1 = round(ep1) + 1;
ep2 = round(ep2) + 1;

% % round the fiber endpoints to align w/ ROIs
% % in = an acpc fiber coordinate
% 
% % round up
% ep_up = floor(ep) + ceil( (ep - floor(ep)) / 0.25) * 0.25;
% 
% % round down
% ep_dn = floor(ep) + floor( (ep - floor(ep)) / 0.25) * 0.25;
% 

%% assign fiber endpoints to labels

% pull fe.roi.coords and convert to acpc coords
fe_roi = mrAnatXformCoords(fe.life.xform.img2acpc, fe.roi.coords);
fe_roi = round(fe_roi) + 1;

% pull all labels (non-zero)
labels = unique(aparc.data);
labels = labels(labels > 0);

display(['Matching streamlines to ' num2str(length(labels)) ' ROIs...']);

% preallocate outputs
out = cell(length(labels), 1);
oout = cell(length(labels), 1);
tfib = zeros(length(labels), 1);
ofib = zeros(length(labels), 1);

% for every label, assign endpoints
tic;
for ii = 1:length(labels)
    
    % catch label info
    out{ii}.label = labels(ii);
    oout{ii}.label = labels(ii);
    
    % pull indices for a label in image space
    [ x, y, z ] = ind2sub(size(aparc.data), find(aparc.data == labels(ii)));
    imgCoords = [ x, y, z ];
    
    % convert label indices to ACPC coordinates
    %acpcCoords = mrAnatXformCoords(fe.life.xform.img2acpc, imgCoords); % rely on fe structure?
    acpcCoords = mrAnatXformCoords(aparc_img2acpc, imgCoords); % rely on file header?
    acpcCoords = round(acpcCoords) + 1;
   
    % catch size 
    out{ii}.size = size(acpcCoords, 1);
    oout{ii}.size = size(acpcCoords, 1);
    
    % test - other out - do I find more streamlines w/o limiting to fe_roi intersection?
    % just find intersections of all fiber end points directly w/ ROIs
    oep1 = ismember(ep1, acpcCoords, 'rows');
    oep2 = ismember(ep2, acpcCoords, 'rows');
    oout{ii}.fibers.end = [ find(oep1); find(oep2) ];
    
    % find label intersection w/ wm mask
    roi_wmi = ismember(fe_roi, acpcCoords, 'rows');
    
    if sum(roi_wmi) <= 0
        warning(['ROI label ' num2str(labels(ii)) ' never intersects the white matter.']);
        out{ii}.fibers.int = [];
        out{ii}.fibers.end = [];
        continue
    end
    
    % subset tensor to label voxels for fiber indices
    [ inds, ~ ] = find(fe.life.M.Phi(:, find(roi_wmi), :));
    
    if size(inds, 2) < 3
        warning(['The found streamline dimensions in label ' num2str(labels(ii)) ' are wrong, none can be identified.']);
        out{ii}.fibers.int = [];
        out{ii}.fibers.end = [];
        continue
    end
    
    % fiber indices of ROI
    out{ii}.fibers.int = unique(inds(:, 3));
    
    % determine if endpoint is in ROI
    
    % subset streamline endpoints to fibers that intersect ROI
    roi_ep1 = ep1(out{ii}.fibers.int, :);
    roi_ep2 = ep2(out{ii}.fibers.int, :);
    
    % for every fiber that intersects, determine if it ends in the label
    frst = ismember(roi_ep1, acpcCoords, 'rows');
    last = ismember(roi_ep2, acpcCoords, 'rows');
            
    % catch and merge the indices of the end points
    out{ii}.fibers.end = [ out{ii}.fibers.int(frst); out{ii}.fibers.int(last) ];
    
    if isempty(out{ii}.fibers.end)
        warning(['No streamlines terminate in ROI ' num2str(labels(ii)) '. Label voxels intersect white matter ROI in ' num2str(sum(roi_wmi)) ' voxels.']);
    end
    
    % total fibers assigned to an endpoint
    tfib(ii) = length(out{ii}.fibers.end);
    ofib(ii) = length(oout{ii}.fibers.end);
    
end
time = toc; 

display(['Successfully assigned endpoints to ' num2str(sum(tfib)) ' of ' num2str(size(fe.life.M.Phi, 3)) ' streamlines in ' num2str(round(time)) ' seconds.']);
display(['Successfully assigned other endpoints to ' num2str(sum(ofib)) ' of ' num2str(size(fe.life.M.Phi, 3)) ' streamlines.']);
clear ii x y z imgCoords acpcCoords roi_wmi inds roi_ep1 roi_ep2 frst last time
clear oep1 oep2

%% build cound of endpoint values

% build every combination
pairs = nchoosek(1:length(labels), 2);

% preallocate paired connection object
pconn = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);
ocon = zeros(length(pairs), 1);

% for every pair of nodes, estimate the connection
for ii = 1:length(pairs)
    
    % create shortcut names
    roi1 = out{pairs(ii, 1)};
    roi2 = out{pairs(ii, 2)};
    oroi1 = oout{pairs(ii, 1)};
    oroi2 = oout{pairs(ii, 2)};
    
    % catch shortcut names
    pconn{ii}.roi1 = roi1.label;
    pconn{ii}.roi2 = roi2.label;
    
    % assign intersections of terminating streamlines
    pconn{ii}.end = intersect(roi1.fibers.end, roi2.fibers.end);
    pconn{ii}.ond = intersect(oroi1.fibers.end, oroi2.fibers.end);
   
    % keep running total of total streamlines assigned a connection
    tcon(ii) = size(pconn{ii}.end, 1);
    ocon(ii) = size(pconn{ii}.ond, 1);
    
end

display(['Build paired connections object with ' num2str(sum(tcon)) ' fibers.']);
display(['Build other paired connections object with ' num2str(sum(ocon)) ' fibers.']);

clear ii roi1 roi2 oroi1 oroi2

%% build matrices from pconn

display('Building Adjacency Matrices...');

% initialize output
emat = zeros(length(labels));
omat = zeros(length(labels));

% for every paired connection
for ii = 1:length(pconn)
   
    % assign pconn values into adjacency matrix
    emat(pairs(ii, 1), pairs(ii, 2)) = size(pconn{ii}.end, 1);
    emat(pairs(ii, 2), pairs(ii, 1)) = size(pconn{ii}.end, 1);
    
    % assign pconn values into other adjacency matrix
    omat(pairs(ii, 1), pairs(ii, 2)) = size(pconn{ii}.ond, 1);
    omat(pairs(ii, 2), pairs(ii, 1)) = size(pconn{ii}.ond, 1);
    
end

clear ii

end

%% check for duplicates

% % for every combination of edges
% comb = nchoosek(1:length(pconn), 2);
% 
% % for every connection
% for ii = 1:length(comb)
%     % the intersection of fiber indices should be empty. if they are not
%     % empty, that's bad. If it's bad, print a warning.
%     bad = ~isempty(intersect(pconn{comb(ii, 1)}.end, pconn{comb(ii, 2)}.end));
%     if bad
%         warning(['Fibers counted in 2 connections. comb index: ' num2str(ii)]);
%     end
% end
% % this should never print a warning if things are working correctly


%% debug plotting

% % make some edge fiber groups
% rfg = fgCreate('name', 'tmp', 'fibers', {fe.fg.fibers{rcon{150}.end}}');
% ifg = fgCreate('name', 'tmp', 'fibers', {fe.fg.fibers{icon{150}.end}}');

% excellent...
%tmpfg = fgCreate('name', 'lrg', 'fibers', {fe.fg.fibers{pconn{962}.end}}');
%bsc_quickPlot(tmpfg, [ 0 0 1 ]);

% % quick plots
% bsc_quickPlot(rfg, [ 1 0 0 ]);
% bsc_quickPlot(ifg, [ 0 0 1 ]);
