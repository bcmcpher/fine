function [ emat, pconn, out ] = feConnectomeSummary(fe, aparc)
%feConnectomeSummary create adjacency matrices from a fit fe structure and a
%    parcellation of cortical regions. 
%   
%   TODO:
%   - split into multiple functions - create pconn and add to that anything I need 
%   - add cleaning
%   - add tract profiles if cell array of labels / files provided
%
% feFile = 'test/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat';
% rois = 'test/inflated_labels.nii.gz';
%

%% load data

% if string is passed, assume it's a path and load it
if isstring(fe) 
    display('Loading fe data...');
    load(fe);
end

% if string is passed, assume it's a path and load it
if isstring(aparc)
    aparc = niftiRead(aparc);
end

%% extract fibers to acpc space and identify endpoint coordinates

display('Coverting streamlines and ROIs to ACPC space...')

% catch xform matrices for aparc
aparc_img2acpc = aparc.qto_xyz;
%aparc_acpc2img = aparc.qto_ijk;

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

% round all the endpoints
ep1 = round(ep1) + 1;
ep2 = round(ep2) + 1;

%% assign fiber endpoints to labels

% % pull fe.roi.coords and convert to acpc coords - do I even need this?
% fe_roi = mrAnatXformCoords(fe.life.xform.img2acpc, fe.roi.coords);
% fe_roi = round(fe_roi) + 1;

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
    acpcCoords = mrAnatXformCoords(aparc_img2acpc, imgCoords); % rely on file header
    acpcCoords = round(acpcCoords) + 1;
   
    % catch size 
    out{ii}.size = size(unique(acpcCoords, 'rows'), 1);
    
    % find streamline endpoints in ROI acpc coordinates
    roi_ep1 = ismember(ep1, acpcCoords, 'rows');
    roi_ep2 = ismember(ep2, acpcCoords, 'rows');
    
    % for fibers that end in rois, catch indices / lengths / weights 
    fibers = [ find(roi_ep1); find(roi_ep2) ];
    out{ii}.end.fibers = fibers;
    out{ii}.end.length = fibLength(out{ii}.end.fibers);
    out{ii}.end.weight = fe.life.fit.weights(out{ii}.end.fibers);
    
    % create endpoint density ROI object
    % create ROI centroid points for glass brain
   
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
%         out{ii}.feroi.int = [];
%         out{ii}.feroi.end = [];
%         continue
%     end
%     
%     % label roi indicies
%     roi_ind = find(roi_wmi);
%     
%     % subset tensor to label voxels for fiber indices
%     [ inds, ~ ] = find(fe.life.M.Phi(:, roi_ind, :));
%     
%     % if inds is not a sufficiently sized input below, go to next ROI
%     if size(inds, 2) < 3
%         warning(['The found streamline dimensions in label ' num2str(labels(ii)) ' are wrong, none can be identified.']);
%         out{ii}.feroi.int = [];
%         out{ii}.feroi.end = [];
%         continue
%     end
%     
%     % fiber indices of ROI
%     out{ii}.feroi.int = unique(inds(:, 3));
%     
%     % determine if endpoint is in ROI
%     
%     % subset streamline endpoints to fibers that intersect ROI
%     feroi_ep1 = ep1(out{ii}.feroi.int, :);
%     feroi_ep2 = ep2(out{ii}.feroi.int, :);
%     
%     % for every fiber that intersects, determine if it ends in the label
%     frst = ismember(feroi_ep1, acpcCoords, 'rows');
%     last = ismember(feroi_ep2, acpcCoords, 'rows');
%             
%     % catch and merge the indices of the end points
%     out{ii}.feroi.end = [ out{ii}.feroi.int(frst); out{ii}.feroi.int(last) ];
%     
%     if isempty(out{ii}.feroi.end)
%         warning(['No streamlines terminate in ROI ' num2str(labels(ii)) '. Label voxels intersect white matter ROI in ' num2str(sum(roi_wmi)) ' voxels.']);
%     end
    
    % total fibers assigned to an endpoint
    tfib(ii) = length(out{ii}.end.fibers);
    
end
time = toc; 

display(['Successfully assigned ' num2str(sum(tfib)) ' of ' num2str(2*size(fe.life.M.Phi, 3)) ' endpoints in ' num2str(round(time)) ' seconds.']);

clear ii x y z imgCoords acpcCoords roi_ep1 roi_ep2 fibers tfib time

%% build  paired connections object

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
    
    % catch size of ROIs
    pconn{ii}.roi1sz = roi1.size;
    pconn{ii}.roi2sz = roi2.size;
    
    % assign intersections of terminating streamlines
    pconn{ii}.all.indices = intersect(roi1.end.fibers, roi2.end.fibers);
    pconn{ii}.all.lengths = fibLength(pconn{ii}.all.indices);
    pconn{ii}.all.weights = fe.life.fit.weights(pconn{ii}.all.indices);
    
    % build matrix weights here
    
    % keep running total of total streamlines assigned a connection
    tcon(ii) = size(pconn{ii}.all.indices, 1);
    
end

display(['Build paired connections object with ' num2str(sum(tcon)) ' fibers.']);

clear ii roi1 roi2 tcon

%% create cleaned version of each connection

% preallocate cleaned edge structure
cln = cell(length(pconn), 1);
clncnt = 0;

% cleaning arguments - set up as arguments?
%maxVolDist = 3;
minLength = 10;
maxDist = 3;
maxLengthStd = 3;
numNodes = 100;
maxIter = 10;

tic;
parfor ii = 1:length(pconn)
    
    if ~isempty(pconn{ii}.all.indices)
        
        % pull indices of connection
        cln{ii}.indices = pconn{ii}.all.indices;
        cln{ii}.lengths = pconn{ii}.all.lengths;
        cln{ii}.weights = pconn{ii}.all.weights;
        
        % drop indices less than a particular length
        cln{ii}.lenindx = cln{ii}.lengths > minLength;
        cln{ii}.indices = cln{ii}.indices(cln{ii}.lenindx);
        cln{ii}.lengths = cln{ii}.lengths(cln{ii}.lenindx);
        cln{ii}.weights = cln{ii}.weights(cln{ii}.lenindx);
        
        % create an fg group of the connection
        cln{ii}.connfib = fgExtract(fg, cln{ii}.indices, 'keep');
        
        % for all streamlines, compute outliers
        [ cln{ii}.afg, cln{ii}.all ] = mbaComputeFibersOutliers(cln{ii}.connfib, maxDist, maxLengthStd, numNodes, 'mean', 0, maxIter);

        % catch image coords for link network
        %cln{ii}.aroi = fefgGet(cln{ii}.afg, 'unique image coords');
        
        % create an fg group for the non-zero weighted fibers
        cln{ii}.nzindx = cln{ii}.weights > 0;
        cln{ii}.nzwfibs = cln{ii}.indices(cln{ii}.nzindx);
        cln{ii}.nzwleng = cln{ii}.lengths(cln{ii}.nzindx);
        cln{ii}.nzwwght = cln{ii}.weights(cln{ii}.nzindx);
        
        % if all streamlines are removed, skip
        if (sum(cln{ii}.nzwfibs > 0) == 0)
            cln{ii}.nzwfibs = double.empty();
            cln{ii}.nzwleng = double.empty();
            cln{ii}.nzwwght = double.empty();
            cln{ii}.nzcnfib = {};
            cln{ii}.nzw = [];
            continue
        end
        
        % pull non-zero streamlines and remove outliers
        cln{ii}.nzcnfib = fgExtract(fg, cln{ii}.nzwfibs, 'keep');
        [ cln{ii}.zfg, cln{ii}.nzw ] = mbaComputeFibersOutliers(cln{ii}.nzcnfib, maxDist, maxLengthStd, numNodes, 'mean', 0, maxIter);

        % catch image coords for link network
        %cln{ii}.zroi = fefgGet(cln{ii}.zfg, 'unique image coords');
        
%         % create connection fg structure
%         cln{ii}.connfib = fgExtract(fg, cln{ii}.indices, 'keep');
%         
%         % pull path neighborhood connection and ROI
%         cln{ii}.pnfibs = feGet(fe, 'Path Neighborhood', cln{ii}.indices);
%         cln{ii}.indroi = feGet(fe, 'coords from fibers', cln{ii}.pnfibs);
%         
%         % initial cleaning
%         cln{ii}.clip = feClipFibersToVolume(cln{ii}.pnfibs, cln{ii}.indroi, maxVolDist);
%         cln{ii}.clip.fibers = mbaFiberSplitLoops(cln{ii}.clip.fibers);
%         
%         % remove fibers below a certain length
%         % how can I keep indices tracked through here?
%         
%         % compute outliers
%         [ ~, keep ] = mbaComputeFibersOutliers(tafg, 3, 3, 100, 'mean', 0, 10);

        % keep count of cleaned connections
        clncnt = clncnt + 1;
        
    else 
        % fill in empty cells and skip to next connection
        cln{ii}.indices = [];
        cln{ii}.lengths = [];
        cln{ii}.weights = [];
        cln{ii}.nzwfibs = [];
        cln{ii}.nzwleng = [];
        cln{ii}.nzwwght = [];
        continue
    end
    
end
time = toc;

display(['Cleaned outlier fibers from ' num2str(clncnt) ' edges in ' num2str(round(time)/60) ' minutes.']);

% add cleaned streamlines to pconn
for jj = 1:length(pconn)
    
    % catch all cleaned indices
    pconn{jj}.cln.indices = cln{jj}.indices;
    pconn{jj}.cln.lengths = cln{jj}.lengths;
    pconn{jj}.cln.weights = cln{jj}.weights;

    % catch non-zero cleaned indices
    pconn{jj}.cln.nzwfibs = cln{jj}.nzwfibs;
    pconn{jj}.cln.nzwleng = cln{jj}.nzwleng;
    pconn{jj}.cln.nzwwght = cln{jj}.nzwwght;
    
end
    
clear ii jj cln clncnt time

%% run parallelized virtural lesion

% preallocate virtual lesion output
vlout = cell(length(pconn), 1);
vlcnt = 0;

tic;
parfor ii = 1:length(pconn)
    
    if sum(pconn{ii}.all.weights > 0) == 0
        
        % set to zero and continue
        vlout{ii}.s.mean  = 0;
        vlout{ii}.em.mean = 0;
        vlout{ii}.j.mean  = 0;
        vlout{ii}.kl.mean = 0;
        continue
                
    else
        % compute a virtual lesion
        [ ewVL, ewoVL ] = feComputeVirtualLesion(fe, pconn{ii}.indices);
        vlout{ii} = feComputeEvidence(ewoVL, ewVL);
        vlcnt = vlcnt + 1;

    end
    
    % add cleaned connections?
    
end
time = toc;

display(['Computed ' num2str(vlcnt) ' vitual lesions in ' num2str(round(time)/60) ' minutes.']);

% add virtual lesion output to pconn
for ii = 1:length(pconn)
    pconn{ii}.all.vl = vlout{ii};
    pconn{ii}.all.matrix.s  = vlout{ii}.s.mean;
    pconn{ii}.all.matrix.em = vlout{ii}.em.mean;
    pconn{ii}.all.matrix.kl = mean(vlout{ii}.kl.mean);
    pconn{ii}.all.matrix.j  = mean(vlout{ii}.j.mean);
end

clear ii vlout vlcnt time

%% build tract profiles

% for ii = 1:length(ms_files)
%     msprop = niftiRead(ms_files{ii});
%     parfor jj = 1:length(pconn)

%% build matrices from pconn

display('Building Adjacency Matrices...');

% initialize output
emat = zeros(length(labels), length(labels), 12);

% for every paired connection
for ii = 1:length(pconn)
    
    % values that can be reused - count, combined voxel size, average length
    cnt = size(pconn{ii}.all.indices, 1);
    psz = pconn{ii}.roi1sz + pconn{ii}.roi2sz;
    len = mean(pconn{ii}.all.lengths);
    
    % redo the traditional count measures w/ non-zero weighted streamlines
    nzw = pconn{ii}.all.weights > 0;
    nzcnt = size(pconn{ii}.all.indices(nzw), 1);
    nzlen = mean(pconn{ii}.all.lengths(nzw));
    
    % add cleaned connections
    
    % build adjacency matrices
    
    % 1. count of streamlines
    emat(pairs(ii, 1), pairs(ii, 2), 1) = cnt;
    emat(pairs(ii, 2), pairs(ii, 1), 1) = cnt;
    pconn{ii}.all.matrix.count = cnt;
    
    % 2. density of streamlines
    dns = (2 * cnt) / psz;
    emat(pairs(ii, 1), pairs(ii, 2), 2) = dns;
    emat(pairs(ii, 2), pairs(ii, 1), 2) = dns;
    pconn{ii}.all.matrix.density = dns;
    
    % 3. length of streamlines
    emat(pairs(ii, 1), pairs(ii, 2), 3) = len;
    emat(pairs(ii, 2), pairs(ii, 1), 3) = len;
    pconn{ii}.all.matrix.length = len;
    
    % 4. density controlling for length (Hagmann, 2008)
    dln = (2 / psz) * sum(1 / pconn{ii}.lengths);
    emat(pairs(ii, 1), pairs(ii, 2), 4) = dln;
    emat(pairs(ii, 2), pairs(ii, 1), 4) = dln;
    pconn{ii}.all.matrix.dnleng = dln;
    
    % 5. non-zero count
    emat(pairs(ii, 1), pairs(ii, 2), 5) = nzcnt;
    emat(pairs(ii, 2), pairs(ii, 1), 5) = nzcnt;
    pconn{ii}.all.matrix.nzcnt = nzcnt;
    
    % 6. non-zero density
    nzdns = (2 * nzcnt) / psz;
    emat(pairs(ii, 1), pairs(ii, 2), 6) = nzdns;
    emat(pairs(ii, 2), pairs(ii, 1), 6) = nzdns;
    pconn{ii}.all.matrix.nzdns = nzdns;
    
    % 7. non-zero length
    emat(pairs(ii, 1), pairs(ii, 2), 7) = nzlen;
    emat(pairs(ii, 2), pairs(ii, 1), 7) = nzlen;
    pconn{ii}.all.matrix.nzlen = nzlen;
    
    % 8. non-zero denstiy * length
    if isempty(pconn{ii}.lengths(nzw)) % if there are no nz lengths
        nzdln = 0;
    else
        nzdln = (2 / psz) * sum(1 / pconn{ii}.lengths(nzw));
    end
    emat(pairs(ii, 1), pairs(ii, 2), 8) = nzdln;
    emat(pairs(ii, 2), pairs(ii, 1), 8) = nzdln;
    pconn{ii}.all.matrix.nzdln = nzdln;
    
    % 9. strength of evidence
    emat(pairs(ii, 1), pairs(ii, 2), 9) = pconn{ii}.matrix.s;
    emat(pairs(ii, 2), pairs(ii, 1), 9) = pconn{ii}.matrix.s;
    
    % 10. strength of evidence
    emat(pairs(ii, 1), pairs(ii, 2), 10) = pconn{ii}.matrix.em;
    emat(pairs(ii, 2), pairs(ii, 1), 10) = pconn{ii}.matrix.em;
    
    % 11. kullback-liebler
    emat(pairs(ii, 1), pairs(ii, 2), 11) = pconn{ii}.matrix.kl;
    emat(pairs(ii, 2), pairs(ii, 1), 11) = pconn{ii}.matrix.kl;
    
    % 12. jefferies divergence
    emat(pairs(ii, 1), pairs(ii, 2), 12) = pconn{ii}.matrix.j;
    emat(pairs(ii, 2), pairs(ii, 1), 12) = pconn{ii}.matrix.j;
    
end

clear ii cnt psz len nzw nzcnt nzlen dns dln nzdns nzdln

% fix impossible values
emat(isnan(emat)) = 0;
emat(isinf(emat)) = 0;
emat(emat < 0) = 0;

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

% excellent...
%tmpfg = fgCreate('name', '', 'fibers', {fe.fg.fibers{pconn{962}.end}}');
%bsc_quickPlot(tmpfg, [ 0 0 1 ]);

%% how many endpoints get found in all ROI labels - became enpoint heatmat plot

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
