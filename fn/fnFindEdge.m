function [ ofg, hmap1, hmap2, conn ] = fnFindEdge(netw, fg, idx1, idx2, smooth)
%fnFindConnection() returns a fiber group between to desired indices if one
% exists.
%
% add profile, matrix, etc. outputs
%

% turn off smoothing if it's not requested
if(~exist('smooth', 'var') || isempty(smooth))
    smooth = [];
end

% IF IDX1 / IDX2 ARE STRINGS...

% sort ROI indices
rois = sort([ idx1 idx2 ]);

% find the edge index of the ROIs
[ ~, idx ] = intersect(netw.parc.pairs, rois, 'rows');
 
% pull the edge
conn = netw.edges{idx};

% inform and exit if empty
if isempty(conn.fibers.indices)
    warning('The requested connection is empty.');
    ofg = [];
    hmap1 = [];
    hmap2 = [];
    conn = [];
    return
    
else
    
    % pull streamline count
    nfib = size(conn.fibers.indices, 1);
    
    % pull roi centers in acpc space
    roi1 = netw.nodes{rois(1)}.centroid.acpc';
    roi2 = netw.nodes{rois(2)}.centroid.acpc';
    
    % create the output connection and length
    ofg = fgCreate('name', [ 'fg_' num2str(rois(1)) '_' num2str(rois(2)) ], ...
                   'colorRgb', [ .1 .25 .65 ], 'fibers', fg.fibers(conn.fibers.indices));
    
    % convert fg into voxel coordinates
    vfg = dtiXformFiberCoords(ofg, netw.parc.xform.acpc2img, 'img');
    
    % reorient fiber group for consistent endpoints - resample b/c I have to
    [ sfg, ep1 ] = dtiReorientFibers(vfg, 100);
    
    % find the distance between the start of the profile and each roi center
    ep_roi1 = norm(ep1 - roi1);
    ep_roi2 = norm(ep1 - roi2);
    
    % if roi2 is closer to the start of the tract profile than roi1, lrflip the fibers
    if (ep_roi2 < ep_roi1)
        ofg.fibers = cellfun(@(x) fliplr(x), ofg.fibers, 'UniformOutput', false);
    end
    
    % preallocate endpoints
    epi = nan(nfib, 3);
    epj = nan(nfib, 3);
    
    % grab all the endpoints in separate matrices
    for fib = 1:nfib
        epi(fib,:) = round(sfg.fibers{fib}(:,1)') + 1;
        epj(fib,:) = round(sfg.fibers{fib}(:,end)') + 1;
    end
    
    % create empty data space
    img1 = zeros(netw.parc.dsize);
    img2 = zeros(netw.parc.dsize);
    
    % create a density map of endpoints for i/j rois
    for ii = 1:nfib
        img1(epi(ii,1), epi(ii, 2), epi(ii, 3)) = img1(epi(ii, 1), epi(ii, 2), epi(ii, 3)) + 1;
        img2(epj(ii,1), epj(ii, 2), epj(ii, 3)) = img2(epj(ii, 1), epj(ii, 2), epj(ii, 3)) + 1;
    end
    
    % smooth data if requested
    if ~isempty(smooth)
        disp(['Smoothing output with a [ ' num2str(smooth) ' ] gaussian kernel...']);
        img1 = smooth3(img1, 'gaussian', smooth);
        img2 = smooth3(img2, 'gaussian', smooth);
    end
    
    % create output niftis in fe space
    hmap1 = niftiCreate('data', img1, 'fname', [ 'ep1_' num2str(rois(1)) '_' num2str(rois(2)) ], ...
                        'qto_xyz', netw.parc.xform.img2acpc, 'qto_ijk', netw.parc.xform.acpc2img);
    
    hmap2 = niftiCreate('data', img2, 'fname', [ 'ep2_' num2str(rois(1)) '_' num2str(rois(2)) ], ...
                        'qto_xyz', netw.parc.xform.img2acpc, 'qto_ijk', netw.parc.xform.acpc2img);
    
    % print some information
    disp(['The edge between ROI''s numbered ' num2str(rois(1)) ' and ' num2str(rois(2)) ' was found at index: ' num2str(idx) ]);
    disp(['This fascicle has ' num2str(size(conn.fibers.indices, 1)) ' streamlines with an average length of ' num2str(mean(conn.fibers.lengths)) ' mm.']);
    
end

end