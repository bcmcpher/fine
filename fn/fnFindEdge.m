function [ ofg, hmap1, hmap2, conn ] = fnFindEdge(netw, fg, idx1, idx2, smooth)
%fnFindEdge() returns a fiber group between to desired indices if one
% exists.
%
% add profile, matrix, etc. outputs
%

%% parse arguments

% turn off smoothing if it's not requested
if(~exist('smooth', 'var') || isempty(smooth))
    smooth = [];
end

%% deal with numeric indices vs. string labels

% pull all labels numbers and names
alabs = cellfun(@(x) x.label, netw.nodes);
aname = cellfun(@(x) x.name, netw.nodes);

% check that names are actually present in netw nodes
notLabeled = any(cellfun(@isempty, aname));

% if they're character inputs and any labels are not present
if (ischar(idx1) || ischar(idx2)) && notLabeled
    error('The input netw is missing node labels, so edges cannot be searched by labels.');
end

% parse the first input
if ischar(idx1) % if the input is character
    [ ~, x ] = intersect(aname, idx1);
else % otherwise look through labels
    [ ~, x ] = intersect(alabs, idx1);
end
if isempty(x) % if it's not found, error
    error('The input ''%d'' is not found in netw.nodes{}', idx1);
end
    
% parse the second input
if ischar(idx2) % if the input is character
    [ ~, y ] = intersect(aname, idx2);
else % otherwise look though labels
    [ ~, y ] = intersect(alabs, idx2);

end
if isempty(y) % if it's not found, error
    error('The input ''%d'' is not found in netw.nodes{}', idx2);
end

%% recover the index into netw.parc.pairs

% sort ROI indices
rois = sort([ x y ]);

% find the edge index of the ROIs
[ ~, idx ] = intersect(netw.parc.pairs, rois, 'rows');
 
%% pull the edge and if it's not empty, return fg / endpoints / edge

conn = netw.edges{idx};

% inform and exit if empty
if isempty(conn.fibers.indices)
    disp(['The edge between ROI''s indexed ' num2str(rois(1)) ' and ' num2str(rois(2)) ' was found at index: ' num2str(idx) ]);
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
    roi1 = netw.nodes{rois(1)}.center.acpc';
    roi2 = netw.nodes{rois(2)}.center.acpc';
    
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
    
    % create usefully labeled nifti objects is labels exist
    if notLabeled
        fname1 = [ 'ep1_' num2str(netw.nodes{rois(1)}.label) ];
        fname2 = [ 'ep2_' num2str(netw.nodes{rois(2)}.label) ];
    else
        fname1 = [ 'ep1_' num2str(netw.nodes{rois(1)}.label) '-' strrep(netw.nodes{rois(1)}.name, ' ', '') ];
        fname2 = [ 'ep2_' num2str(netw.nodes{rois(2)}.label) '-' strrep(netw.nodes{rois(2)}.name, ' ', '') ];
    end
    
    % create output niftis in fe space
    hmap1 = niftiCreate('data', img1, 'fname', fname1, ...
                        'qto_xyz', netw.parc.xform.img2acpc, 'qto_ijk', netw.parc.xform.acpc2img);
    
    hmap2 = niftiCreate('data', img2, 'fname', fname2, ...
                        'qto_xyz', netw.parc.xform.img2acpc, 'qto_ijk', netw.parc.xform.acpc2img);
    
    % print some information
    disp(['The edge between ROI''s indexed ' num2str(rois(1)) ' and ' num2str(rois(2)) ' was found at index: ' num2str(idx) ]);
    disp(['This fascicle has ' num2str(size(conn.fibers.indices, 1)) ' streamlines with an average length of ' num2str(mean(conn.fibers.lengths)) ' mm.']);
    
end

end