function [ netw ] = fnAveragePropertyEdges(netw, fg, vol, dtype, keep, clobv, clobd)
%fnAverageEdgeProperty finds the central tendency of a tract from an 
% appropriately aligned volume or dt6 structure. Does not require ENCODE
%
% INPUTS:
%     netw  - is the network object wanting edge volumes and central tendency.
%
%     fg    - the whole brain tractography in ac-pc space
%
%     vol   - a volume of data or a dt6 on which to store data and compute
%             central tendency per edge. Must be in alignment with fe ROI 
%             for valid results.
%
%     dtype - a field name the central tendency can be stored with
%
%     keep  - whether or not to store coordinates / values for each edge
%
%     clobv - whether or not to overwrite volume data if the fields already exist 
%
%     clobd - whether or not to overwrite dtype central tendency if the fields already exist 
%
% OUTPUTS:
%     netw  - the network object with the voxel indices, volume, and central tendency
%             for each edge computed based on the mapping of fg to vol
%
% NOTES:
% - allow a list of input volumes / labels? - check dimensions to do multiple at once?
% - coords for wmvol are available but not stored - it's big and probably has no use.
% - coords for edges are available but not stored, used for link networks.
% - - add it back when I get to those? Or just compute as needed?
% - vals for dtype are available but not stored, used in MI comparison?
% - - add back when they're necessary, useful to have but maybe not needed.
% - store more specific central tendency (standard error, etc.)
%
% EXAMPLE:
%
% % load data
% parc          = niftiRead('labels.nii.gz');
% favol         = niftiRead('fa.nii.gz');
% fg            = feGet(fe, 'fibers acpc');
% fibers        = fg.fibers;
% fibLength     = fefgGet(fg, 'length');
% weights       = feGet(fe, 'fiberweights');
% Phi           = feGet(fe, 'Phi');
% roi           = feGet(fe, 'roicoords');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % find unique voxels in each edge and add them to pconn
% pconn = fnFindPathVoxels(pconn, 'nzw', Phi);
%
% % compute central tendency of external volume
% pconn = fnAverageEdgeProperty(pconn, 'nzw', roi, favol);
%
% Brent McPherson (c), 2017 - Indiana University
%

if(~exist('keep', 'var') || isempty(keep))
    keep = false;
end

if(~exist('clobv', 'var') || isempty(clobv))
    clobv = 0;
end

if(~exist('clobd', 'var') || isempty(clobd))
    clobd = 0;
end

if(isfield(vol, 'dt6'))
    
    disp('Computing edge central tendency from dt6 structure...');
    
    tol = 1e-6;
    dt = vol;
    
    % extract tensor data from dt6
    [ ~, eigVal ] = dtiEig(dt.dt6);
    mask = repmat(all(eigVal == 0, 4),[ 1 1 1 3 ]);
    eigVal(eigVal < tol) = tol;
    eigVal(mask) = 0;
    clear mask;
    
    switch dtype
        
        case 'fa'
            dat = dtiComputeFA(eigVal);
            
        case 'md'
            [ ~, dat ] = dtiComputeFA(eigVal);
            
        case 'rd'
            [ ~, ~, dat ] = dtiComputeFA(eigVal);
                    
        case 'ad'
            [ ~, ~, ~, dat ] = dtiComputeFA(eigVal);
    
        case 'cl'
            dat = dtiComputeWestinShapes(eigVal, 'lsum');
            
        case 'cp'
            [ ~, dat ] = dtiComputeWestinShapes(eigVal, 'lsum');
            
        case 'cs'
            [ ~, ~, dat ] = dtiComputeWestinShapes(eigVal, 'lsum');
            
        otherwise
            error('A volume summary of ''%s'' cannot be extracted from a dt6', dtype);
    end
    
    % create the "nifti" volume of the data in AC-PC space
    vol = niftiCreate('data', dat, ...
                      'qto_xyz', dt.xformToAcpc);
                  
    % create vol filename
    name = [ 'dt6_' dtype ];

else
    disp('Computing edge central tendency from volume...');
    
    % create vol filename
    [ ~, nam, ext ] = fileparts(vol.fname);
    name = [ nam ext ];

end

% extract relavant data from volume
data = vol.data;

% pull xforms to move b/w fibers acpc / micro image space
micro_acpc2img = niftiGet(vol, 'qto_ijk');

% grab volume voxel size for scaling volume
volvx = prod(vol.pixdim);

% grab the volume size
volsz = size(vol.data);

%% check the size of the image passed to see if it matches parcellation

% print the longer warning
be_careful = false;

% clearly state if this volume is sliced different than the parcellation
if ~all(netw.parc.dsize == volsz)
    warning('The image dimensions of the volume do not match the current dimensions.');
    be_careful = true;
end

% clearly state if this volume has a different voxel size than the parcellation
if ~all(prod(netw.parc.voxmm) == volvx)
    warning('The voxel dimensions of the volume do not match the current dimensions.');
    be_careful = true;
end

% If any of those are true, there is no guarantee that the volume passed will
% reliably map to all the stored central tendency values. The streamlines
% can be mapped to different resolutions, and as long as they have a small
% enough step size they should be mapped and estimated correctly. However,
% this will result in the estimated volume being different between the 
% dtypes, which makes some sanity checks fail. It is good practice to only
% use data in the same resolution as the diffusion data for parcellations
% and edge property estimates.
if be_careful
    warning('The volume that will be sampled is not sliced the same as the rest of the data. If this assumption is violated, the edge volumes may not map exactly to all stored dtypes.');
end

%% safety checks - exit before loop if an overwrite is expected

% don't overwrite volume if it's already written before clobber check
write_vol = true;

% logic to stop volume rewrite
if isfield(netw.volume, 'wmvol')
    if (isfield(netw.volume, 'wmvol') && clobv ~= 1)
        warning('Volumes for these edges are already computed. Set clobv = 1 to recompute.');
        write_vol = false;
    else      
        warning('Network volume has been previously estimated. clobv set to overwrite previously estimated volume.');   
    end
end

% more strictly stop dtype label redos - this doesn't change on accident
if isfield(netw.edges{1}, 'micros')
    if isfield(netw.edges{1}.micros, dtype)
        if clobd ~= 1
            error('Central tendencies for this edge label are already computed. Set clobd = 1 to recompute.');
        else
            warning('Edge central tendency has been previously estimated. clobd set to overwrite previously estimated values.');
        end
    end
end

%% create whole brain summary

fprintf('Estimating white matter volume...\n');

% transform whole brain fg into voxel space
vfg = dtiXformFiberCoords(fg, micro_acpc2img, 'img');

% pull the total number of nodes
MAXMEM = 4000000; % max of 4gb?
mnodes = round(MAXMEM*15000000/32000000);
tnodes = size(horzcat(vfg.fibers{:}), 2);
nfiber = size(vfg.fibers, 1);

% split into batches like LiFE encoding
nBatch = max([5,ceil(tnodes/mnodes)]); % number of batches
mFiber = ceil(nfiber/nBatch); % number of fibers per batch

% preallocate array for unique voxels in each streamline
ufg = cell(nfiber, 1);

% for every batch
fprintf('Converting streamlines to voxels in batches...\n');
for batch = 1:nBatch
    fprintf('Encoding batch %02.f\n',batch)
    
    % pull the fiber range to not explode memory
    fibers_range = (batch-1)*mFiber + 1: min(batch*mFiber,nfiber);
    
    % cell array to convert img space streamline nodes to the unique voxel indices
    bat = cellfun(@(x) unique((round(x)+1)', 'rows'), vfg.fibers(fibers_range), 'UniformOutput', false);
    
    % cell array to fix converted voxel indices to be within volume indices
    out = cellfun(@fixVoxelBounds, bat, repmat(mat2cell(volsz, 1, 3), [ size(fibers_range, 2) 1 ]), 'UniformOutput', false); 
    % fixVoxelBounds defined below
    
    % catch the output for later use
    ufg(fibers_range) = out;
end

% find all unique voxels across all streamlines
coords = unique(vertcat(ufg{:}), 'rows');

% % replace 0 indices w/ 1 - sanity check
% coords(coords <= 1) = 1;
% 
% % replace indices exceeding max w/ max index - sanity check
% coords(coords(:, 1) > volsz(1), 1) = volsz(1);
% coords(coords(:, 2) > volsz(2), 2) = volsz(2);
% coords(coords(:, 3) > volsz(3), 3) = volsz(3);
% % these fixes match how ENCODE checks fg

clear batch bat

% vectorized call through streamlines rounding nodes to voxel indices
%stvx = cellfun(@(x) unique((round(x)+1)', 'rows'), vfg.fibers, 'UniformOutput', false);
%coords = unique(vertcat(stvx{:}), 'rows'); 
% THIS SHOULD BE TOO BIG FOR MEMORY - SO CHUNK LIKE LiFE CODE
% this matches the LiFE voxel ROI if the volume passed is sliced like DWI

% estimate white matter volume in structure volume
wmvox = size(coords, 1);
wmvol = wmvox * volvx;
netw.volume.wmvol = wmvol;
if keep % optionally store full wm ROI
    netw.volume.coords = coords;
end
clear coords % it's big otherwise...

fprintf('Total white matter volume estimated at: %d mm^3.\n', wmvol);
fprintf('Finding volume and central tendency of ''%s'' for all edges...\n', name);

tic;
for ii = 1:length(netw.edges)
    
    % pull the connection
    conn = netw.edges{ii};
    
    % if the connection is empty, fill in zeros
    if isempty(conn.fibers.indices)
        
        % fill in empty voxel coords
        conn.volume.volume = 0;
        conn.volume.prpvol = 0;
        if keep
            conn.micros.(dtype).raw = [];
        end
        conn.micros.(dtype).mn = 0;
        conn.micros.(dtype).md = 0;
        conn.micros.(dtype).sd = 0;
        
    else
        
        % otherwise pull the voxels for the tract to make a volume
        tcoords = unique(vertcat(ufg{conn.fibers.indices}), 'rows');
        
        % optionally store the voxel indices for each tract
        if keep
            conn.volume.coords = tcoords;
        end
        
        % extract the micro data for the edge
        vals = nan(size(tcoords, 1), 1);
        for coord = 1:size(tcoords, 1)
            vals(coord) = data(tcoords(coord, 1), tcoords(coord, 2), tcoords(coord, 3));
        end
        
        % store the volume
        if write_vol
            conn.volume.volume = size(vals, 1) * volvx;
            conn.volume.prpvol = conn.volume.volume / wmvol;
        end

        % optionally store raw values and central tendency measures
        if keep
            conn.micros.(dtype).raw = vals;
        end
        conn.micros.(dtype).mn = mean(vals, 'omitnan');
        conn.micros.(dtype).md = median(vals, 'omitnan');
        conn.micros.(dtype).sd = std(vals, 'omitnan');
                
    end
    
    % store connection back in output
    netw.edges{ii} = conn;
        
end
time = toc;

disp([ 'Found the volume and central tendency of ' name ' for all edges in ' num2str(time) ' seconds.' ]);

end

% a function to vectorize the correction of rounded voxels to be w/in volume indices
function [ coords ] = fixVoxelBounds(coords, volsz)

% replace 0 indices w/ 1 - sanity check
coords(coords <= 1) = 1;

% replace indices exceeding max w/ max index - sanity check
coords(coords(:, 1) > volsz(1), 1) = volsz(1);
coords(coords(:, 2) > volsz(2), 2) = volsz(2);
coords(coords(:, 3) > volsz(3), 3) = volsz(3);
% these fixes match how ENCODE checks fg

end