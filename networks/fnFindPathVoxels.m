function [ pconn ] = fnFindPathVoxels(Phi, pconn, label)
%fnFindPathVoxels finds the tensor indices of voxels for every connection 
% to create either microstructural summaries of networks or for
% generating link networks.
%
% INPUTS:
%     Phi   - the sparse tensor object from LiFE model to find voxel /
%             streamline node intersection.
%
%     pconn - is the paired connections object to compute edge volumes for.
%
%     label - string indicating the fiber groups for which to create virtual lesions
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%             Additionally, this can be run after cleaning, resulting in
%             valid calls of 'all_clean' and 'nzw_clean', respectively.
%
% OUTPUTS:
%     omat - 3d array containing (nodes x nodes x edge_type) of processed networks
%     olab - cell array of labels of the 'edge_type' along the 3rd dimension of omat
%
% TODO:
% - compute proportion of white matter (unique_edge / voxel dim of Phi)
%
% EXAMPLE:
%
% % load data
% parc          = niftiRead('labels.nii.gz');
% fg            = feGet(fe, 'fibers acpc');
% fibers        = fg.fibers;
% fibLength     = fefgGet(fg, 'length');
% weights       = feGet(fe, 'fiberweights');
% Phi           = feGet(fe, 'Phi');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % find unique voxels in each edge and add them to pconn
% pconn = fnFindPathVoxels(Phi, pconn, 'nzw');
%
% Brent McPherson (c), 2017 - Indiana University
%

display('Finding voxels in tracts for all edges...')

tic;
parfor ii = 1:length(pconn)
    
    % pull the connection
    conn = pconn{ii}.(label);
    
    % if the connection is empty, fill in zero
    if isempty(conn.indices)
        
        % fill in empty voxel coords
        conn.pvoxels = [];
    
    else
        
        % pull subtensor of the connection
        subtensor = Phi(:, :, conn.indices);
        
        % pull the unique voxel indices of the connection
        conn.pvoxels = unique(subtensor.subs(:, 2));
    
    end
    
    % re-assign connection
    pconn{ii}.(label) = conn;
    
end
time = toc;

display([ 'Found the voxels for all edges in ' num2str(time/60) ' minutes.' ]);

end

