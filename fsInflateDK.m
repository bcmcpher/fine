function [ fh, ilabs, wm ] = fsInflateDK(aparc, infl, nlog)
%% 
% currently assigns mode value of neighbors
% can repeat inflation, unsure of ideal optimum
% what else should I compute?
% decide what files to save - prompts / dummy checks
% ROI centroids for glass brain
% check to see if I can import FS directly?
%
% aparc = '105115/mri/aparc+aseg.nii.gz'; infl = 3;

%% import and set up data

display('Loading aparc+aseg nifti...');

% read in aligned aparc+aseg .nii (read labels)
rlabs = niftiRead(aparc);

%% create brain mask

mask = rlabs;
mask.fname = 'brain_mask.nii.gz';
mask.data = mask.data > 0;

%% create wm mask from aparc+aseg labels - from Franco

wm = rlabs;
wm.fname = 'wm_mask.nii.gz';

% labels for wm mask to expand into
invals = [2 41 16 17 28 60 51 53 12 52 13 18 54 50 11 251 252 253 254 255 10 49 46 7];

% dummy counters
wmCounter = 0; noWMCounter = 0;

% create the wm mask
origvals = unique(wm.data(:));
for ii = 1:length(origvals);
    if any(origvals(ii) == invals)
        wm.data( wm.data == origvals(ii) ) = 1;
        wmCounter=wmCounter+1;
    else
        wm.data( wm.data == origvals(ii) ) = 0;
        noWMCounter = noWMCounter + 1;
    end
end

display(['White matter masked with ' num2str(wmCounter) ' regions; ' num2str(noWMCounter) ' regions remain.']);

%% grab labels for inflation

% subset to original labels
olabs = rlabs;
olabs.fname = 'original_labels.nii.gz';

% set labels less than 1000 to 0
olabs.data(olabs.data < 1000) = 0;

% remove unknown labels
olabs.data(olabs.data == 1000 | olabs.data == 2000) = 0;

% compute centroids of ROIs
% find unique labels / isolate voxels, averaged x/y/z to return

%% probably drop labels mask...

% make orig label mask
mlabs = olabs;
mlabs.fname = 'orig_labels_mask.nii.gz';
mlabs.data(mlabs.data > 0) = 1;

%% create local neighborhood logicals

% can probably streamline these better...

% vertices are all points around voxel, 26 neigborhood
ivert = ones(3, 3, 3);

% edges share an border boundary, 18 neighborhood
iedge = ivert;
iedge(1, 1, 1) = 0; iedge(1, 3, 1) = 0; iedge(3, 1, 1) = 0; iedge(3, 3, 1) = 0;
iedge(1, 1, 3) = 0; iedge(1, 3, 3) = 0; iedge(3, 1, 3) = 0; iedge(3, 3, 3) = 0;

% faces share boundary with voxel, 6 neighborhood
iface = iedge;
iface(2, 1, 1) = 0; iface(1, 2, 1) = 0; iface(3, 2, 1) = 0; iface(2, 3, 1) = 0;
iface(2, 1, 3) = 0; iface(1, 2, 3) = 0; iface(3, 2, 3) = 0; iface(2, 3, 3) = 0;
iface(1, 1, 2) = 0; iface(1, 3, 2) = 0; iface(3, 1, 2) = 0; iface(3, 3, 2) = 0;

switch nlog
    case{'vert'}
        nmsk = ivert;
    case{'edge'}
        nmsk = iedge;
    case{'face'}
        nmsk = iface;
    otherwise
        nmsk = ivert;
        warning('Invalid neighborhood logical requested. Can either select: ''vert'', ''edge'', or ''face''. Defaulting to ''vert''');
end

%% pad data set and loop over every entry

display(['Inflating ' num2str(length(unique(olabs.data(:)))-1) ' ROIs...']);

% catch data in empty array of original size
out = olabs.data;

% padded original dataset for first iteration
pdat = padarray(olabs.data, [1, 1, 1]);

% loop for number of inflations
for hh = 1:infl
    
    display(['Running ROI inflation pass: ' num2str(hh)]);
    
    % for every voxel in the data
    for ii = 1:size(olabs.data, 1)
        for jj = 1:size(olabs.data, 2)
            for kk = 1:size(olabs.data, 3)
                
                % if outside brain mask, skip
                if mask.data(ii, jj, kk) == 0
                    continue; % if it's in the mask...
                end
                
                % create padded indices
                pii = ii + 1;
                pjj = jj + 1;
                pkk = kk + 1;
                
                % if plabs is labeled skip b/c I don't change labels
                if pdat(pii, pjj, pkk) > 0 % and if it's not labeled...
                    continue;
                end
                                
                % find neighborhood
                neigh = pdat(pii-1:pii+1, pjj-1:pjj+1, pkk-1:pkk+1);
                
                % logical index of face / edge / vert
                neigh = neigh(nmsk == 1);
                
                % only keep neighbors > 0
                neigh = neigh(neigh > 0);
                
                % figure out the most popular neighbor
                lab = mode(neigh);
                
                % assign label based on mode if voxel is 0 and not assigning 0
                if out(ii, jj, kk) == 0 && lab ~= 0
                    
                    % output label is assigned
                    out(ii, jj, kk) = lab;
                    
                end
            end
        end
    end
    
    % recompute search space after every iteration
    pdat = padarray(out, [1, 1, 1]);
    
end

% create inflated labels nifti
ilabs = olabs;
ilabs.fname = 'inflated_labels.nii.gz';
ilabs.data = out;

display(['Saving inflated ROIs in: ' ilabs.fname '...']);

% save inflated labels
niftiWrite(ilabs, ilabs.fname);

%% debug figure

% pull a slice and see what is added
x = olabs.data(:,:,78);
y = ilabs.data(:,:,78);
z = y - x;
a = (ilabs.data(:,:,78) > 0) + (mask.data(:,:,78) > 0);
wmOverlap1 = wm.data(:,:,78);
wmOverlap2 = z > 0;
wmOverlap = wmOverlap1 + wmOverlap2;

fh = figure; 
subplot(3, 2, 1); imagesc(x); title('Original Labels'); subplot(3, 2, 2); imagesc(y); title('Inflated Labels'); 
subplot(3, 2, 3); imagesc(a); title('Labels + Mask'); subplot(3, 2, 4); imagesc(mask.data(:,:,78)>0); title('Mask');
subplot(3, 2, 5); imagesc(z); title('Label Expansion'); subplot(3, 2, 6); imagesc(wmOverlap); title('WM Intersect');

end