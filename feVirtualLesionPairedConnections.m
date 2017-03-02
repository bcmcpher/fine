function [ pconn ] = feVirtualLesionPairedConnections(pconn, label)
%feVirtualLesionPairedConnections runs virtual lesion on a field of pconn indices 
%   
% feFile = 'test/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat';
% rois = 'test/inflated_labels.nii.gz';
%


% run parallelized virtural lesion

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
