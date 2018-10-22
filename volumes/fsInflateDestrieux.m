function [ ilabs, wm, dbout ] = fsInflateDestrieux(aparc, infl, nlog, subc, outfile)
%% 
%
% check to see if I can import FS directly?
%
% aparc = '105115/mri/aparc.a2009s+aseg.nii.gz'; infl = 3;
%

%% parse optional argument for including (or not) subcortical labels

if(~exist('subc', 'var') || isempty(subc))
    subc = 0;
end

%% import and set up data

display('Loading aparc.a2009.aseg nifti...');

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
for ii = 1:length(origvals)
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

% predefined cortex labels for this parcellation
ctxvals = [ 11101:11175 12101:12175 ];

% predefined extra 7 subcortical in each hemi
sbcvals = [ 10 11 12 13 17 18 26 49 50 51 52 53 54 58 ];

% brainstem / cerebellum?

% missing labels of useful but otherwise covered (?) regions
% l/r thalamus: 9 48
% l/r insula 19 55
% l/r opperculum 20 56
% l/r substantia nigra: 27 59

% subset to original labels
olabs = rlabs;
olabs.fname = 'original_labels.nii.gz';

% if subcortical labels are requested
if subc == 1
    
    display('Selecting cortical and subcortical labels for dilation...');
    
    % merge all labels
    vals = [ sbcvals ctxvals ];

    % copy logic of white matter mask, but fill in labels
    ovals = unique(olabs.data(:));
    for ii = 1:length(ovals)
        if any(ovals(ii) == vals)
            olabs.data( olabs.data == ovals(ii) ) = ovals(ii);
        else
            olabs.data( olabs.data == ovals(ii) ) = 0;
        end
    end
    
else
    
    display('Selecting cortical labels for dilation...');
    
    % just use cortical labels
    vals = ctxvals; 
    
    % copy logic of white matter mask, but fill in labels
    ovals = unique(olabs.data(:));
    for ii = 1:length(ovals)
        if any(ovals(ii) == vals)
            olabs.data( olabs.data == ovals(ii) ) = ovals(ii);
        else
            olabs.data( olabs.data == ovals(ii) ) = 0;
        end
    end

    % set labels less than 1000 to 0
    olabs.data(olabs.data < 1000) = 0;
    
end

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
        prop = 7;
    case{'edge'}
        nmsk = iedge;
        prop = 6;
    case{'face'}
        nmsk = iface;
        prop = 2;
    otherwise
        nmsk = ivert;
        warning('Invalid neighborhood logical requested. Can either select: ''vert'', ''edge'', or ''face''. Defaulting to ''vert''');
end

%% pad data set and loop over every entry

display(['Inflating ' num2str(length(unique(olabs.data(:)))-1) ' ROIs...']);

% catch data in empty array of original size
out = olabs.data;

%% create debug counter for ROIs to determine when they stop growing

% process to get count of each label
dbcnt = sort(out(out > 0));
p = find([numel(dbcnt); diff(dbcnt); numel(dbcnt)]);
label = dbcnt(p(1:end-1))';
count = diff(p)';

% intitialize returned counter
dbout = zeros(infl + 1, length(label));
dbout(1, :) = label;
dbout(2, :) = count;

%% resume fxn

% padded original dataset for first iteration
pdat = padarray(out, [1, 1, 1]);

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
                if pdat(pii, pjj, pkk) > 0 % if it's not labeled...
                    continue;
                end
                                
                % find neighborhood
                neigh = pdat(pii-1:pii+1, pjj-1:pjj+1, pkk-1:pkk+1);
                
                % logical index of face / edge / vert
                neigh = neigh(nmsk == 1);
                
                % only keep of prop # of neighbors > 0
                
                if sum(neigh > 0) > prop
                    lab = mode(neigh(neigh > 0));
                else
                    lab = 0;
                end
                
                % don't assign zero
                if lab == 0 
                    continue;
                end
                
                % assign label based on mode if voxel is 0 is still white matter
                if out(ii, jj, kk) == 0 && wm.data(ii, jj, kk) == 1
                    
                    % output label is assigned
                    out(ii, jj, kk) = lab;
                    
                end
            end
        end
    end
    
    % recompute search space after every iteration
    pdat = padarray(out, [1, 1, 1]);
    
    % process to get count of each label
    dbcnt = sort(out(out > 0));
    p = find([numel(dbcnt); diff(dbcnt); numel(dbcnt)]);
    count = diff(p)';
    
    % update debug out count
    dbout(hh+2,:) = count;
    
end

% create inflated labels nifti
ilabs = olabs;
ilabs.fname = outfile;
ilabs.data = out;

display(['Saving inflated ROIs in: ' ilabs.fname '...']);

% save inflated labels
niftiWrite(ilabs, ilabs.fname);

end
