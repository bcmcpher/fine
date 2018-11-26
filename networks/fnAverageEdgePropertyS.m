function [ pconn ] = fnAverageEdgePropertyS(pconn, label, fg, vol, dtype, meas, clobber)
%fnAverageEdgePropertyS finds the central tendency of a tract from an 
% appropriately aligned volume or dt6 structure. Does not require ENCODE
%
% INPUTS:
%     pconn - is the paired connections object to compute edge volumes for.
%
%     label - string indicating the fiber groups for which to create virtual lesions
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%             Additionally, this can be run after cleaning, resulting in
%             valid calls of 'all_clean' and 'nzw_clean', respectively.
%
%     fg    - the streamlines in ac-pc space
%
%     vol   - a volume of data or a dt6 on which to store data and compute 'meas'
%             per edge. Must be in alignment with fe ROI for valid results
%
%     dtype - a field name the central tendency can be stored with
%
%     meas  - the measure of central tendency that will summarize data
%             options are:
%                 - mean
%                 - median
%                 - std
%                 - var
%
%     clobber - whether or not to overwrite data if the field(s) already exist 
%
% OUTPUTS:
%     pconn - the paired connection object with the voxel indices for each 
%             edge computed from the fe.roi
%
% TODO:
% - if computing a central tendency on stored data, not require 
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

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

if(isfield(vol, 'dt6'))
    
    disp('Computing tract central tendency from dt6...');
    
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
    disp('Computing tract central tendency from volume...');
    
    % create vol filename
    [ ~, nam, ext ] = fileparts(vol.fname);
    name = [ nam ext ];

end

% extract relavant data from volume
data = vol.data;
    
% output field name
field = [ meas '_' dtype ];

% pull xforms to move b/w fibers acpc / micro image space
micro_acpc2img = niftiGet(vol, 'qto_ijk');

fprintf('Finding ''%s'' of ''%s'' for all edges...\n', meas, name);

tic;
for ii = 1:length(pconn)
    
    % pull the connection
    conn = pconn{ii}.(label);
    
    % if the connection is empty, fill in zero
    if isempty(conn.indices)
        
        % fill in empty voxel coords
        conn.(dtype).raw = [];
        conn.(dtype).(meas) = 0;
        vals = [];
        
    else
        
        % if raw is already stored and clobber is blank, use stored values
        if (isfield(conn, [ 'raw_' dtype ]) && clobber == 0)
            
            %vals = conn.([ 'raw_' dtype ]);
            vals = conn.(dtype).raw;
            
        else
            
            % create the fascicle subset
            tfg = fg;
            tfg.fibers = fg.fibers(conn.indices);
            
            % convert to image space
            tfg = dtiXformFiberCoords(tfg, micro_acpc2img, 'img');
            
            % create one array of streamline nodes
            snd = horzcat(tfg.fibers{:})';
            
            % create the mask
            coords = unique(round(snd), 'rows');
            
            % replace 0 indices w/ 1
            coords(coords <= 0) = 1;
            
            % replace indices exceeding max w/ max index
            mnmx = size(vol.data);
            coords(coords(:, 1) > mnmx(1), 1) = mnmx(1);
            coords(coords(:, 2) > mnmx(2), 2) = mnmx(2);
            coords(coords(:, 3) > mnmx(3), 3) = mnmx(3);
            
            % dev - code to dilate volume for central tendency
            % makes edges much larger (~4x) and MUCH, MUCH slower (1 min to 4 hours)
            % would this be the volume of edges as well?
            
%             % create a mask volume of tract
%             tvol = zeros(mnmx);
%             for jj = 1:size(coords, 1)
%                 tvol(coords(jj, 1), coords(jj, 2), coords(jj, 3)) = 1; 
%             end
%             % dilate / fill / erode mask for a more even volume?
%             se1 = strel('cube', 2);
%             se2 = strel('cube', 1);
%             ivol = imdilate(tvol, se1); 
%             ivol = imfill(ivol, 'holes');
%             ivol = imerode(ivol, se2);
%             % pull vals for volume for central tendency
%             vals = vol.data(ivol == 1);
            
            % extract the micro data for the edge
            vals = nan(size(coords, 1), 1);
            for coord = 1:size(coords, 1)
                vals(coord) = data(coords(coord, 1), coords(coord, 2), coords(coord, 3));
            end

            % add the raw values
            conn.(dtype).raw = vals;
            
        end
        
        % compute the central tendency
        switch meas
            case {'mean', 'average'}
                conn.(dtype).(meas) = mean(vals, 'omitnan');
            case {'median'}
                conn.(dtype).(meas) = median(vals, 'omitnan');
            case {'std', 'sd'}
                conn.(dtype).(meas) = std(vals, 'omitnan');
            case {'var', 'variance'}
                conn.(dtype).(meas) = var(vals, 'omitnan');
            otherwise
                error('''%s'' central tendency is not currently defined for this context.', meas);
        end

    end
    
    % re-assign connection and matrix value
    pconn{ii}.(label).(dtype).vals = vals;
    pconn{ii}.(label).(dtype).(meas) = conn.(dtype).(meas);
    pconn{ii}.(label).matrix.(field) = conn.(dtype).(meas);
    
end
time = toc;

disp([ 'Found the ' meas ' of ' name ' for all edges in ' num2str(time) ' seconds.' ]);

end
