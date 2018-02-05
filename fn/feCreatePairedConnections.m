function [ pconn, rois, out ] = feCreatePairedConnections(parc, fibers, fibLength, weights)
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

%% extract fibers to acpc space and identify endpoint coordinates

display('Assigning streamline endpoints to ROI labels...')

% catch xform matrices for parc
parc_acpc2img = niftiGet(parc, 'qto_ijk');
parc_img2acpc = niftiGet(parc, 'qto_xyz');

% convert acpc fibers to parcellation space
fg = fgCreate;
fg.fibers = fibers;
fg = dtiXformFiberCoords(fg, parc_acpc2img, 'img');

% re-extract label space fiber coordinates
fibers = fg.fibers;

% initialize endpoint outputs
ep1 = zeros(length(fibers), 3);
ep2 = zeros(length(fibers), 3);

% for every fiber, pull the end points
if size(fibers{1}, 1) ~= 3
    error('Expected fibers with size(3, N)'); 
end

% extract every streamline end point
parfor ii = 1:length(fibers)
    ep1(ii,:) = fibers{ii}(:,1)';
    ep2(ii,:) = fibers{ii}(:,end)';
end

% compute size of fg for reporting
fibers_size = size(fibers, 1);

clear ii fibers

% round converted end points to match image indices
% produces highest assignment compared to floor / ceil / round() + 1
ep1 = round(ep1);
ep2 = round(ep2);

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
parfor ii = 1:length(labels)
    
    % catch label info
    rois{ii}.label = labels(ii);
    
    % pull indices for a label in image space
    [ x, y, z ] = ind2sub(parcsz, find(parc_data == labels(ii)));
    imgCoords   = [ x, y, z ];

    % catch size of ROI
    rois{ii}.size = size(unique(imgCoords, 'rows'), 1);
    
    % find streamline endpoints in image coordinates for label
    roi_ep1 = ismember(ep1, imgCoords, 'rows');
    roi_ep2 = ismember(ep2, imgCoords, 'rows');
    
    % convert end points in ROI to streamline indices
    fibers = [ find(roi_ep1); find(roi_ep2) ];
    
    % for fibers that end in rois, catch indices / lengths / weights
    rois{ii}.end.fibers = fibers;
    rois{ii}.end.length = fibLength(rois{ii}.end.fibers);
    rois{ii}.end.weight = weights(rois{ii}.end.fibers);
    
    % create endpoint ROI in vistasoft format
    rois{ii}.roi = dtiNewRoi(num2str(labels(ii)), 'red', [ ep1(roi_ep1, :); ep2(roi_ep2, :) ]);
    
    % create ROI centroid
    rois{ii}.centroid.img = round(mean(imgCoords) + 1);
    rois{ii}.centroid.acpc = mrAnatXformCoords(parc_img2acpc, rois{ii}.centroid.img);
    
    if isempty(rois{ii}.end.fibers)
        warning(['ROI label ' num2str(labels(ii)) ' has no streamline terminations.']);
    end
    
    % total fibers assigned to an endpoint
    tfib(ii) = length(rois{ii}.end.fibers);
    
end

display(['Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*fibers_size) ' endpoints.']);

clear ii x y z imgCoords roi_ep1 roi_ep2 fibers time

%% build paired connections object

% build every unique combination of labels
pairs = nchoosek(1:length(labels), 2);

% preallocate paired connection object
pconn = cell(length(pairs), 1);
tcon = zeros(length(pairs), 1);
ncon = zeros(length(pairs), 1);

display('Building paired connections...');

% for every pair of nodes, estimate the connection
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
    
    % find all weighted streamlines
    nzw = pconn{ii}.all.weights > 0;
    
    pconn{ii}.nzw.indices = pconn{ii}.all.indices(nzw);
    pconn{ii}.nzw.lengths = pconn{ii}.all.lengths(nzw);
    pconn{ii}.nzw.weights = pconn{ii}.all.weights(nzw);

    % find remaining zero weighted streamlines (for debugging)
    zrw = pconn{ii}.all.weights == 0;

    pconn{ii}.zrw.indices = pconn{ii}.all.indices(zrw);
    pconn{ii}.zrw.lengths = pconn{ii}.all.lengths(zrw);
    pconn{ii}.zrw.weights = pconn{ii}.all.weights(zrw);

    % calculate combined size of ROI
    psz = pconn{ii}.roi1sz + pconn{ii}.roi2sz;
    
    % values used to calculate different edge weights
    cnt = size(pconn{ii}.all.indices, 1);
    len = nanmean(pconn{ii}.all.lengths);
    dln = sum(1 / pconn{ii}.all.lengths);
    
    % values used to calculate different edge weights for non-zero weighted streamlines
    nzcnt = size(pconn{ii}.all.indices(nzw), 1);
    nzlen = nanmean(pconn{ii}.all.lengths(nzw));
    
    % if there are no nz lengths
    if isempty(pconn{ii}.nzw.lengths) 
        nzdln = 0;
    else
        nzdln = sum(1 ./ pconn{ii}.nzw.lengths);
    end
    
    % create equivalent measures for zero weighted fibers
    zrcnt = size(pconn{ii}.all.indices(zrw), 1);
    zrlen = nanmean(pconn{ii}.all.lengths(zrw));
    zrdln = sum(1 ./ pconn{ii}.zrw.lengths);
    
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
    
    % create zero weighted streamline counts
    pconn{ii}.zrw.matrix.count = zrcnt;
    pconn{ii}.zrw.matrix.density = (2 * zrcnt) / psz;
    pconn{ii}.zrw.matrix.length = zrlen;
    pconn{ii}.zrw.matrix.denlen = (2 / psz) * zrcnt * zrdln;
        
    % keep running total of total streamlines assigned a connection
    tcon(ii) = size(pconn{ii}.all.indices, 1);
    
    % if the connection is empty or not
    if size(pconn{ii}.all.indices, 1) > 0
        ncon(ii) = 1;
    end
    
end

display(['Built paired connections object with ' num2str(sum(tcon)) ' streamlines.']);

%% simple summary of counts

% total streamlines, # of end points assigned
out.strmct = fibers_size;
out.assign = sum(tfib);

% total connections, # assigned connections, # of streamlines assinged to connections
out.nconns = size(pconn, 1);
out.aconns = sum(ncon);
out.paired = sum(tcon);

% parcellation voxel size
out.dvoxmm = abs(parc_img2acpc([1, 6, 11]));

end
