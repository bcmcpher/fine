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
    clobv = false;
end

if(~exist('clobd', 'var') || isempty(clobd))
    clobd = false;
end

% assume only 1 input
ninputs = 1;

% pull the number of edges
nedges = length(netw.edges);

% check if there are cell array of input names / volumes
if iscell(dtype) || iscell(vol)
    % check if they're the same length
    if (length(dtype) == length(vol))
        % if they labels / volumes have the same length, check dimensions
        dims = cellfun(@(x) size(x.data), vol, 'UniformOutput', false);
        if isequal(dims{:})
            disp('Array of inputs have matching dimensions.');
        else
            warning('Dimensions of input volumes do not match. This may cause issues.');
        end
        ninputs = length(dtype);
    else
        error('The cell array of labels and volumes are not the same length. Unable to appropriately assign inputs.');
    end
else
    
    % make the single inputs cells for indexing sanity
    dtype = {dtype};
    vol = {vol};
    
end

% check for dt6, but not if an array of inputs is passed
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
            error('A volume summary of ''%s'' cannot be extracted from a dt6', dtype{1});
    end
    
    % create the "nifti" volume of the data in AC-PC space
    vol = {niftiCreate('data', dat, 'qto_xyz', dt.xformToAcpc)};
                  
    % create vol filename
    name = {[ 'dt6_' dtype{1} ]};

else
%     disp('Computing edge central tendency from volume...');
%     
%     % create vol filename
%     [ ~, nam, ext ] = fileparts(vol.fname);
%     name = [ nam ext ];

end

% preallocate inputs
name = cell(ninputs, 1);
data = cell(ninputs, 1);
micro_acpc2img = cell(ninputs, 1);
volvx = cell(ninputs, 1);
volsz = cell(ninputs, 1);

% for every input, pull data
for input = 1:ninputs
    
    % create vol filename
    [ ~, nam, ext ] = fileparts(vol{input}.fname);
    name{input} = [ nam ext ];
    
    % extract relavant data from volume
    data{input} = vol{input}.data;
    
    % pull xforms to move b/w fibers acpc / micro image space
    micro_acpc2img{input} = niftiGet(vol{input}, 'qto_ijk');
    
    % grab volume voxel size for scaling volume
    volvx{input} = prod(vol{input}.pixdim);
    
    % grab the volume size
    volsz{input} = size(vol{input}.data);
    
end

% round to "functionally" identical spaces (they should actually be...)
% tensor / noddi are not identical, but are at 4 decimal places. close enough?
micro_acpc2img = cellfun(@(x) round(x, 4), micro_acpc2img, 'UniformOutput', false);

% hard exit if the inputs are not in the same space.
if (ninputs > 1) && ~isequal(micro_acpc2img{:})
    error('The array of inputs are in different spaces. There is no sane way to do this.');
end

% get the number of fibers passed
nfib = size(fg.fibers, 1);

%% check the size of the image passed to see if it matches parcellation

% print the longer warning
be_careful = false;

% clearly state if this volume is sliced different than the parcellation
if ~all(netw.parc.dsize == volsz{1})
    warning('The image dimensions of the new volume(s) do not match the parcellations stored dimensions.');
    be_careful = true;
end

% clearly state if this volume has a different voxel size than the parcellation
if ~all(prod(netw.parc.voxmm) == volvx{1})
    warning('The voxel dimensions of the new volume(s) do not match the parcellations stored dimensions.');
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
    warning('The new volume to be sampled is not sliced the same as the currenlty stored parcellation. If this assumption is violated, edge volumes may not map exactly to all stored dtypes.');
end

%% safety checks - exit before loop if an overwrite is expected

% don't overwrite volume if it's already written before clobber check
write_vol = true;

% logic to stop volume rewrite
if isfield(netw.volume, 'wmvol')
    if (isfield(netw.volume, 'wmvol') && ~clobv)
        warning('Volumes for these edges are already computed. Set clobv = 1 to recompute.');
        write_vol = false;
    else      
        warning('Network volume has been previously estimated. clobv set to overwrite previously estimated volume.');   
    end
end

% more strictly stop dtype label redos - this doesn't change on accident
if isfield(netw.edges{1}, 'micros')
    if isfield(netw.edges{1}.micros, dtype)
        if ~clobd
            error('Central tendencies for this edge label are already computed. Set clobd = 1 to recompute.');
        else
            warning('Edge central tendency has been previously estimated. clobd set to overwrite previously estimated values.');
        end
    end
end

%% create whole brain summary

fprintf('Estimating white matter volume from fg...\n');

% transform whole brain fg into voxel space
% able to use the first index b/c they all have to match or there's only 1
vfg = dtiXformFiberCoords(fg, micro_acpc2img{1}, 'img');

% vectorized call through streamlines rounding nodes to voxel indices
ufg = cellfun(@fixVoxelBounds, vfg.fibers, repmat(mat2cell(volsz{1}, 1, 3), [ nfib 1 ]), 'UniformOutput', false);
coords = unique(vertcat(ufg{:}), 'rows'); 
% THIS COULD BE TOO BIG FOR MEMORY - SO CHUNK LIKE LiFE CODE BELOW
% coords exactly matches the LiFE voxel ROI if the volume passed is sliced like DWI
% this vectorized approach doesn't blow up memory or matlab updated (?)

% % pull the total number of nodes
% MAXMEM = 4000000; % max of 4gb?
% mnodes = round(MAXMEM*15000000/32000000);
% tnodes = size(horzcat(vfg.fibers{:}), 2);
% nfiber = size(vfg.fibers, 1);
% 
% % split into batches like LiFE encoding
% nBatch = max([5,ceil(tnodes/mnodes)]); % number of batches
% mFiber = ceil(nfiber/nBatch); % number of fibers per batch
% 
% % preallocate array for unique voxels in each streamline
% ufg = cell(nfiber, 1);
% 
% % for every batch
% fprintf('Converting streamlines to voxels in batches...\n');
% for batch = 1:nBatch
%     fprintf('Encoding batch %02.f\n',batch)
%     
%     % pull the fiber range to not explode memory
%     fibers_range = (batch-1)*mFiber + 1: min(batch*mFiber,nfiber);
%     
%     % cell array to fix converted voxel indices to be within volume indices
%     out = cellfun(@fixVoxelBounds, vfg.fibers(fibers_range), repmat(mat2cell(volsz, 1, 3), [ size(fibers_range, 2) 1 ]), 'UniformOutput', false); 
%     % fixVoxelBounds defined below
%     
%     % catch the output for later use
%     ufg(fibers_range) = out;
%     
% end
% 
% % find all unique voxels across all streamlines
% coords = unique(vertcat(ufg{:}), 'rows');
% 
% clear batch bat

% estimate white matter volume in structure volume
wmvox = size(coords, 1);
wmvol = wmvox * volvx{1};
netw.volume.wmvol = wmvol;

% optionally store full wm ROI
if keep 
    netw.volume.coords = coords;
end
clear coords % it's too big otherwise...

fprintf('Total white matter volume estimated at: %d mm^3.\n', wmvol);

tic;

fprintf('Computing volume and central tendency within all edges of:\n');
fprintf(' - %s\n', name{:});

for edge = 1:nedges
    
    % pull the connection
    conn = netw.edges{edge};
    
    % if the connection is empty, fill in zeros
    if isempty(conn.fibers.indices)
        
        warning('Edge index %d is empty. This should be impossible.', edge);
        
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
        
        % for every input
        for input = 1:ninputs
        
            % preallocate the empty values
            vals = nan(size(tcoords, 1), 1);
            
            % store the input coords for this edge
            for coord = 1:size(tcoords, 1)
                vals(coord) = data{input}(tcoords(coord, 1), tcoords(coord, 2), tcoords(coord, 3));
            end
        
            % store the volume
            if write_vol && (input == 1) % only do this the first time
                conn.volume.volume = size(vals, 1) * volvx{input};
                conn.volume.prpvol = conn.volume.volume / wmvol;
            end
        
            % optionally store raw values and central tendency measures
            if keep
                conn.micros.(dtype{input}).raw = vals;
            end
            
            % write the central tendency
            conn.micros.(dtype{input}).mn = mean(vals, 'omitnan');
            conn.micros.(dtype{input}).md = median(vals, 'omitnan');
            conn.micros.(dtype{input}).sd = std(vals, 'omitnan');
            
        end
        
    end
    
    % store connection back in output
    netw.edges{edge} = conn;
    
end

time = toc;

disp([ 'Found the volume and central tendency for all edges in ' num2str(time) ' seconds.' ]);

end

%% function to facilitate vectorized call across streamlines

% a function to vectorize the correction of rounded voxels to be w/in volume indices
function [ coords ] = fixVoxelBounds(coords, volsz)

% convert streamline coords in img space to voxel indices like LiFE does
coords = unique((round(coords)+1)', 'rows');

% replace 0 indices w/ 1 - sanity check
coords(coords <= 1) = 1;

% replace indices exceeding max w/ max index - sanity check
coords(coords(:, 1) > volsz(1), 1) = volsz(1);
coords(coords(:, 2) > volsz(2), 2) = volsz(2);
coords(coords(:, 3) > volsz(3), 3) = volsz(3);
% these fixes match how ENCODE checks fg

end
