function [ emat, pconn, out ] = feConnectomeSummary(feFile, rois)
%feConnectomeSummary create adjacency matrices from a fit fe structure and a
%    parcellation of cortical regions. 
%   
%   TODO:
%   - too many streamlines are unassigned
%   - add computations for other matrices
%   - add cleaning
%   - parallelize virtual lesion
%   - add tract profiles if cell array of labels / files provided
%
% feFile = 'test/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat';
% rois = 'test/inflated_labels.nii.gz';
%

%% load data

display('Loading data...');

% load the arguments - build checks to use structures if provided
load(feFile);
aparc = niftiRead(rois);

% catch xform matrices for aparc
aparc_img2acpc = aparc.qto_xyz;
%aparc_acpc2img = aparc.qto_ijk;

%% extract fibers to acpc space and identify endpoint coordinates

% convert fibers to acpc space
fg = feGet(fe, 'fg acpc');

% get fiber lengths
fibLength = fefgGet(fg, 'length');

% % cleaning can start here, but I need with and without - move to
% %     individual connections w/ other steps later...
% % cleaning argument - minimum fiber length
% minLength = 10;
% 
% % go ahead and filter by length - cleaning
% lindx = fibLength > minLength;
% 
% % subset minimum length filtered fibers
% fg.fibers = fg.fibers(lindx);
% fe.life.fit.weights = fe.life.fit.weights(lindx);
% fe.life.M.Phi = fe.life.M.Phi(:,:,find(lindx));
%
% the rest of the cleaning is in plot_path_neighborhood

% initialize endpoint outputs
ep1 = zeros(length(fg.fibers), 3);
ep2 = zeros(length(fg.fibers), 3);

% for every fiber, pull the end points
for ii = 1:length(fg.fibers)
    ep1(ii,:) = fg.fibers{ii}(:,1)';
    ep2(ii,:) = fg.fibers{ii}(:,end)';
end

% round all the endpoints
ep1 = round(ep1) + 1;
ep2 = round(ep2) + 1;

% round the fiber endpoints to align w/ ROIs exactly - bad idea, useful code
% 
% % round up
% ep1 = floor(ep1) + ceil( (ep1 - floor(ep1)) / 0.25) * 0.25;
% ep2 = floor(ep2) + ceil( (ep2 - floor(ep2)) / 0.25) * 0.25;
% 
% % round down
% ep1 = floor(ep1) + floor( (ep1 - floor(ep1)) / 0.25) * 0.25;
% ep2 = floor(ep2) + floor( (ep2 - floor(ep2)) / 0.25) * 0.25;

%% assign fiber endpoints to labels

% pull fe.roi.coords and convert to acpc coords - do I even need this?
%fe_roi = mrAnatXformCoords(fe.life.xform.img2acpc, fe.roi.coords);
%fe_roi = round(fe_roi) + 1;

% pull all labels (non-zero)
labels = unique(aparc.data);
labels = labels(labels > 0);

display(['Matching streamlines to ' num2str(length(labels)) ' nodes...']);

% preallocate outputs
out = cell(length(labels), 1);
tfib = zeros(length(labels), 1);

% for every label, assign endpoints
tic;
for ii = 1:length(labels)
    
    % catch label info
    out{ii}.label = labels(ii);
    
    % pull indices for a label in image space
    [ x, y, z ] = ind2sub(size(aparc.data), find(aparc.data == labels(ii)));
    imgCoords = [ x, y, z ];
    
    % convert label indices to ACPC coordinates
    %acpcCoords = mrAnatXformCoords(fe.life.xform.img2acpc, imgCoords); % rely on fe structure
    acpcCoords = mrAnatXformCoords(aparc_img2acpc, imgCoords); % rely on file header
    acpcCoords = round(acpcCoords) + 1;
   
    % catch size 
    out{ii}.size = size(acpcCoords, 1);
    
    % test - other out - do I find more streamlines w/o limiting to fe_roi intersection?
    % just find intersections of all fiber end points directly w/ ROIs -
    % just as fast, but finds more
    roi_ep1 = ismember(ep1, acpcCoords, 'rows');
    roi_ep2 = ismember(ep2, acpcCoords, 'rows');
    
    % for fibers that end in rois, catch indices / lengths / weights 
    out{ii}.end.fibers = [ find(roi_ep1); find(roi_ep2) ];
    out{ii}.end.length = fibLength(out{ii}.end.fibers);
    out{ii}.end.weight = fe.life.fit.weights(out{ii}.end.fibers);
    
    if isempty(out{ii}.end.fibers)
        warning(['ROI label ' num2str(labels(ii)) ' has no streamline terminations.']);
    end
    
%  %  % this is closer to my original logic - look for intersection in white matter mask
%  %  % this finds fewer, and if wm rois inflate correctly the steps below are irrelevant
%  %  % don't need tensor to align ROIs - unless node order is encoded it's
%  %  % not worth dealing with tensor to assign fibers.
%
%     % find label intersection w/ wm mask
%     roi_wmi = ismember(fe_roi, acpcCoords, 'rows');
%     
%     if sum(roi_wmi) <= 0
%         warning(['ROI label ' num2str(labels(ii)) ' never intersects the white matter.']);
%         out{ii}.fibers.int = [];
%         out{ii}.fibers.end = [];
%         continue
%     end
%     
%     % subset tensor to label voxels for fiber indices
%     [ inds, ~ ] = find(fe.life.M.Phi(:, find(roi_wmi), :));
%     
%     if size(inds, 2) < 3
%         warning(['The found streamline dimensions in label ' num2str(labels(ii)) ' are wrong, none can be identified.']);
%         out{ii}.fibers.int = [];
%         out{ii}.fibers.end = [];
%         continue
%     end
%     
%     % fiber indices of ROI
%     out{ii}.fibers.int = unique(inds(:, 3));
%     
%     % determine if endpoint is in ROI
%     
%     % subset streamline endpoints to fibers that intersect ROI
%     roi_ep1 = ep1(out{ii}.fibers.int, :);
%     roi_ep2 = ep2(out{ii}.fibers.int, :);
%     
%     % for every fiber that intersects, determine if it ends in the label
%     frst = ismember(roi_ep1, acpcCoords, 'rows');
%     last = ismember(roi_ep2, acpcCoords, 'rows');
%             
%     % catch and merge the indices of the end points
%     out{ii}.fibers.end = [ out{ii}.fibers.int(frst); out{ii}.fibers.int(last) ];
%     
%     if isempty(out{ii}.fibers.end)
%         warning(['No streamlines terminate in ROI ' num2str(labels(ii)) '. Label voxels intersect white matter ROI in ' num2str(sum(roi_wmi)) ' voxels.']);
%     end
    
    % total fibers assigned to an endpoint
    tfib(ii) = length(out{ii}.end.fibers);
    
end
time = toc; 

display(['Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*size(fe.life.M.Phi, 3)) ' endpoints in ' num2str(round(time)) ' seconds.']);
clear ii x y z imgCoords acpcCoords roi_ep1 roi_ep2 time

%% build cound of endpoint values

% build every unique combination of labels
pairs = nchoosek(1:length(labels), 2);

% preallocate paired connection object
pconn = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);

% for every pair of nodes, estimate the connection
for ii = 1:length(pairs)
    
    % create shortcut names
    roi1 = out{pairs(ii, 1)};
    roi2 = out{pairs(ii, 2)};
    
    % catch shortcut names
    pconn{ii}.roi1 = roi1.label;
    pconn{ii}.roi2 = roi2.label;
    
    % assign intersections of terminating streamlines
    pconn{ii}.end = intersect(roi1.end.fibers, roi2.end.fibers);
   
    % keep running total of total streamlines assigned a connection
    tcon(ii) = size(pconn{ii}.end, 1);
    
end

display(['Build paired connections object with ' num2str(sum(tcon)) ' fibers.']);

clear ii roi1 roi2 

%% build matrices from pconn

display('Building Adjacency Matrices...');

% initialize output
emat = zeros(length(labels));

% for every paired connection
for ii = 1:length(pconn)
   
    % assign pconn values into adjacency matrix
    emat(pairs(ii, 1), pairs(ii, 2)) = size(pconn{ii}.end, 1);
    emat(pairs(ii, 2), pairs(ii, 1)) = size(pconn{ii}.end, 1);
    
end

clear ii

end

%% how many endpoints get found in all ROI labels

% % find and convert all ROI labels to acpc space
% [ x1, y1, z1 ] = ind2sub(size(aparc.data), find(aparc.data > 0));
% imgLabels = [ x1, y1, z1 ];
% acpcLabel = mrAnatXformCoords(aparc_img2acpc, imgLabels);
% acpcLabel = round(acpcLabel) + 1;
% clear x1 y1 z1 imgLabels
% 
% % total number of first and last endpoints
% sum(ismember(ep1, acpcLabel, 'rows'))
% sum(ismember(ep2, acpcLabel, 'rows'))
% 
% % find indices of end point fibers
% x1 = find(ismember(ep1, acpcLabel, 'rows'));
% y1 = find(ismember(ep2, acpcLabel, 'rows'));
% 
% % find unassinged end point indices
% x2 = find(~ismember(ep1, acpcLabel, 'rows'));
% y2 = find(~ismember(ep2, acpcLabel, 'rows'));
% 
% % find total number of fibers that have both endpoints assigned
% total = intersect(x1, y1);
% 
% % combine assinged / unassigned enpoints
% aep = [ ep1(x1, :); ep2(y1, :) ];
% uep = [ ep1(x2, :); ep2(y2, :) ];
% 
% clear x1 y1 x2 y2
% 
% % pull a random sample of points so it will even plot
% saep = aep(randperm(1000), :);
% suep = uep(randperm(2000), :);
% 
% % quick plot of assigned enpoints
% figure; hold on
% for ii = 1:size(saep, 1)
%     plot3(saep(ii, 1), saep(ii, 2), saep(ii, 3), 'o', 'markersize', 7, 'color', 'blue');
% end
% 
% for ii = 1:size(suep, 1)
%     plot3(suep(ii, 1), suep(ii, 2), suep(ii, 3), '.', 'markersize', 10, 'color', 'red');
% end
% 
% axis equal; axis square;
% view(0, 0);
% clear ii

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

% excellent...
%tmpfg = fgCreate('name', '', 'fibers', {fe.fg.fibers{pconn{962}.end}}');
%bsc_quickPlot(tmpfg, [ 0 0 1 ]);
