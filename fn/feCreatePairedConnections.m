function [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights)
%feCreatePairedConnections creates pconn object of every possible unique pair of 
% of labels in parc w/ streamlines from fe object. 
%   
% TODO:
% - make sure this is as fast as possible without having to run too many extra steps for more basic tasks
% - catch endpoint ROI data - potential cosine similarity comparison between ROIs
%
% feFile = 'test/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat';
% rois = 'test/inflated_labels.nii.gz';
%

%% load data

% % if string is passed, assume it's a path and load it
% if isstring(fe) 
%     display('Loading fe data...');
%     load(fe);
% end

% if string is passed, assume it's a path and load it
if isstring(parc)
    parc = niftiRead(parc);
end

%% extract fibers to acpc space and identify endpoint coordinates

display('Coverting streamlines and ROIs to ACPC space...')

% catch xform matrices for aparc
parc_img2acpc = niftiGet(parc,'qto_xyz');
%parc_acpc2img = niftiGet(parc,'qto_ijk');

% initialize endpoint outputs
ep1 = zeros(length(fibers), 3);
ep2 = zeros(length(fibers), 3);

% for every fiber, pull the end points
% fibers = fg.fibers; clear fg
if size(fibers{1},1)~=3; error('Expected fibers with size(3,N)'); end
parfor ii = 1:length(fibers)
    ep1(ii,:) = fibers{ii}(:,1)';
    ep2(ii,:) = fibers{ii}(:,end)';
end
fibers_size = size(fibers, 1);
clear ii fibers

% Trasform from AC-PC mm coordinates to correspodning voxel (in AC-PC T1w anatomy)
ep1 = round(ep1) + 1;
ep2 = round(ep2) + 1;

%% assign fiber endpoints to labels

% pull all labels (non-zero)
labels = unique(parc.data);
labels = labels(labels > 0);

display(['Matching streamlines to ' num2str(length(labels)) ' nodes...']);

% preallocate outputs
rois = cell(length(labels), 1);
tfib = zeros(length(labels), 1);

parcsz = size(parc.data);
parc_data = parc.data;

% for every label, assign endpoints
tic;

parfor ii = 1:length(labels)
    
    % catch label info
    rois{ii}.label = labels(ii);
    
    % pull indices for a label in image space
    [ x, y, z ] = ind2sub(parcsz, find(parc_data == labels(ii)));
    imgCoords   = [ x, y, z ];
    
    % convert label indices to ACPC coordinates
    acpcCoords = mrAnatXformCoords(parc_img2acpc, imgCoords); % rely on file header
    acpcCoords = round(acpcCoords) + 1;
    
    % catch size
    rois{ii}.size = size(unique(acpcCoords, 'rows'), 1);
    
    % find streamline endpoints in ROI acpc coordinates
    roi_ep1 = ismember(ep1, acpcCoords, 'rows');
    roi_ep2 = ismember(ep2, acpcCoords, 'rows');
    
    % for fibers that end in rois, catch indices / lengths / weights
    fibers = [ find(roi_ep1); find(roi_ep2) ];
    rois{ii}.end.fibers = fibers;
    rois{ii}.end.length = fibLength(rois{ii}.end.fibers);
    rois{ii}.end.weight = weights(rois{ii}.end.fibers);
    
    % create ROI centroid
    rois{ii}.centroid.acpc = round(mean(acpcCoords) + 1);
    rois{ii}.centroid.img  = round(mean(imgCoords)) + 1;
    
%     % create endpoint density ROI object
% 
%     % combine found acpc endpoints
%     all_acpc_ep = [ ep1(roi_ep1, :); ep2(roi_ep2, :) ];
%     
%     % convert acpc points back to image space
%     all_img_ep = round(mrAnatXformCoords(parc_acpc2img, all_acpc_ep));
%     
%     % identify unique voxels and counts in image space
%     [ unq, ~, cnt ] = unique(all_img_ep, 'rows');
%     
%     % convert ACPC end points back to image space voxel coords
%     %unq = round(mrAnatXformCoords(parc_acpc2img, unq));
%     
%     % catch image space coordinates and counts in roi structure 
%     rois{ii}.volume = [ unq accumarray(cnt, 1) ];
    
%     % write these down to check alignment
%     img = zeros(size(parc.data));
%     for jj = 1:size(y{1}.volume, 1)
%         img(y{1}.volume(jj, 1), y{1}.volume(jj, 2), y{1}.volume(jj, 3)) = y{1}.volume(jj, 4); 
%     end
%     z = niftiCreate('data', img, 'fname', 'testROI.nii.gz', 'qto_xyz', parc.qto_xyz, 'qto_ijk', parc.qto_ijk);

    if isempty(rois{ii}.end.fibers)
        warning(['ROI label ' num2str(labels(ii)) ' has no streamline terminations.']);
    end
    
    % total fibers assigned to an endpoint
    tfib(ii) = length(rois{ii}.end.fibers);
    
end
time = toc;

display(['Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*fibers_size) ' endpoints in ' num2str(round(time)/60) ' minutes.']);

clear ii x y z imgCoords acpcCoords roi_ep1 roi_ep2 fibers tfib time

%% build paired connections object

% build every unique combination of labels
pairs = nchoosek(1:length(labels), 2);

% preallocate paired connection object
pconn = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);

display('Building paired connections...');

% for every pair of nodes, estimate the connection
tic;
parfor ii = 1:length(pairs)
    
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
    pconn{ii}.all.weights = weights(pconn{ii}.all.indices);
    
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
    if isempty(  pconn{ii}.nzw.lengths) % if there are no nz lengths
        nzdln = 0;
    else
        nzdln = sum(1 / pconn{ii}.nzw.lengths);
    end
    
    % create all streamline counts
    pconn{ii}.all.matrix.count = cnt;
    pconn{ii}.all.matrix.density = (2 * cnt) / psz;
    pconn{ii}.all.matrix.length = len;
    pconn{ii}.all.matrix.denlen = (2 / psz) * cnt * dln;
        
    % create non-zero weighted streamline counts
    pconn{ii}.nzw.matrix.count = nzcnt;
    pconn{ii}.nzw.matrix.density = (2 * nzcnt) / psz;
    pconn{ii}.nzw.matrix.length = nzlen;
    pconn{ii}.nzw.matrix.denlen = (2 / psz) * nzcnt * nzdln;
    
    % keep running total of total streamlines assigned a connection
    tcon(ii) = size(pconn{ii}.all.indices, 1);
    
end
time = toc;

display(['Built paired connections object with ' num2str(sum(tcon)) ' streamlines in ' num2str(time/60) ' minutes.']);

end
