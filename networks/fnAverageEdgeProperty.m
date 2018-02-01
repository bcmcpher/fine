function [ pconn ] = fnAverageEdgeProperty(pconn, label, roi, vol, dtype, meas, clobber)
%fnAverageEdgeProperty finds the central tendency of a tract from an 
% appropriately aligned volume.
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
%     vol   - a volume of data on which central tendecies will be computed
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
% - second output that holds raw coordinates too?
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

% create vol filename
[ ~, nam, ext ] = fileparts(vol.fname);
name = [ nam ext ];

% output field name
field = [ meas '_' dtype ];

% extract relavant data from volume
data = vol.data;

% error without clobber set if measure is computed
if (isfield(pconn{1}.(label).volume, field) && clobber == 0)
    error('''%s'' of ''%s'' for this label is already computed. Set clobber = 1 to recompute.', meas, dtype);
end

fprintf('Finding ''%s'' of ''%s'' for all edges...\n', meas, name);

tic;
parfor ii = 1:length(pconn)
    
    % pull the connection
    conn = pconn{ii}.(label).volume;
    
    % if the connection is empty, fill in zero
    if isempty(conn.pvoxels)
        
        % fill in empty voxel coords
        conn.([ 'raw_' dtype ]) = [];
        conn.(field) = [];
            
    else
        
        % if raw is already stored and clobber is blank, use stored values
        if (isfield(conn, [ 'raw_' dtype ]) && clobber == 0)
            
            vals = conn.([ 'raw_' dtype ]);
            
        else
            
            % pull coordinates of connection
            coords = roi(conn.pvoxels, :);
            
            % grab the values
            vals = zeros(size(coords, 1), 1);
            
            for jj = 1:size(coords, 1)
                vals(jj) = data(coords(jj, 1), coords(jj, 2), coords(jj, 3));
            end
            
            % add the raw values
            conn.([ 'raw_' dtype ]) = vals;
            
        end
        
        % compute the central tendency
        switch meas
            case {'mean', 'average'}
                conn.(field) = mean(vals);
            case {'median'}
                conn.(field) = median(vals);
            case {'std', 'sd'}
                conn.(field) = std(vals);
            case {'var', 'variance'}
                conn.(field) = var(vals);
            otherwise
                error('''%s'' central tendency is not currently defined for this context.', meas);
        end

    end
    
    % re-assign connection and matrix value
    pconn{ii}.(label).volume = conn;
    pconn{ii}.(label).matrix.(field) = conn.(field);
    
end
time = toc;

display([ 'Found the ' meas ' of ' name ' for all edges in ' num2str(time/60) ' minutes.' ]);

end


