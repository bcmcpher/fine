function [ netw ] = fnAverageEdgeProperty(netw, roi, vol, dtype, meas, clobber)
%%
% REPLACED WITH fnAverageEdgePropertyS() TO NOT REQUIRE ENCODE STRUCTURE
% AND FOLLOW UP WITH fnFindPathVoxels()
%%
%fnAverageEdgeProperty finds the central tendency of a tract from an 
% appropriately aligned volume or dt6 structure.
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
%     roi   - the coordinates of the indices in the volume from the fe
%             structure.
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
%field = [ meas '_' dtype ];

% NEED TO FIX OR REMOVE
% % error without clobber set if measure is computed
% if (isfield(netw.pconn{1}.volume, dtype) && clobber == 0)
%     if (isfield(netw.pconn{1}.volume.(dtype), meas) && clobber == 0)
%         error('''%s'' of ''%s'' is already computed. Set clobber = 1 to recompute.', meas, dtype);
%     end
% end

fprintf('Finding ''%s'' of ''%s'' for all edges...\n', meas, name);

tic;
for ii = 1:length(netw.pconn)
    
    % pull the connection
    conn = netw.pconn{ii};
    
    % if the connection is empty, fill in zero
    if isempty(conn.pvoxels)
        
        % fill in empty voxel coords
        conn.volume.(dtype).raw = [];
        conn.volume.(dtype).(meas) = 0;
        
    else
        
        % if raw is already stored and clobber is blank, use stored values
        if (isfield(conn, [ 'raw_' dtype ]) && clobber == 0)
            
            %vals = conn.([ 'raw_' dtype ]);
            vals = conn.volume.(dtype).raw;
            
        else
            
            % pull coordinates of connection
            coords = roi(conn.pvoxels, :);
            
            % grab the values
            vals = zeros(size(coords, 1), 1);
            
            for jj = 1:size(coords, 1)
                vals(jj) = data(coords(jj, 1), coords(jj, 2), coords(jj, 3));
            end
            
            % debugging - write out the volume that is averaged
            tvol = vol;
            tvol.data = zeros(size(vol.data));   
            for coord = 1:size(coords,1)
                tvol.data(coords(coord, 1), coords(coord, 2), coords(coord, 3)) = vals(coord);
            end
            niftiWrite(tvol, [ '/geode2/home/u010/bcmcpher/Carbonate/fine-data/edges/tedge_' sprintf('%02d',ii) '.nii.gz' ]);
            
            % add the raw values
            conn.volume.(dtype).raw = vals;
            
        end
        
        % compute the central tendency
        switch meas
            case {'mean', 'average'}
                conn.volume.(dtype).(meas) = mean(vals, 'omitnan');
            case {'median'}
                conn.volume.(dtype).(meas) = median(vals, 'omitnan');
            case {'std', 'sd'}
                conn.volume.(dtype).(meas) = std(vals, 'omitnan');
            case {'var', 'variance'}
                conn.volume.(dtype).(meas) = var(vals, 'omitnan');
            otherwise
                error('''%s'' central tendency is not currently defined for this context.', meas);
        end

    end
    
    % re-assign connection and matrix value
    netw.pconn{ii} = conn;
    %pconn{ii}.matrix.(field) = conn.(dtype).(meas);
    
end
time = toc;

disp([ 'Found the ' meas ' of ' name ' for all edges in ' num2str(time) ' seconds.' ]);

end


