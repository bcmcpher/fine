function [ pconn ] = fnFindPathVoxels(pconn, label, Phi, dim, clobber)
%fnFindPathVoxels finds the tensor indices of voxels for every connection 
% to create either microstructural summaries of networks or for
% generating link networks.
%
% INPUTS:
%     pconn   - is the paired connections object to compute edge volumes for.
%
%     label   - string indicating the fiber groups for which to create virtual lesions
%               either:
%                       'all' for all assigned streamlines or
%                       'nzw' for non-zero weighted fibers returned by LiFE
%               Additionally, this can be run after cleaning, resulting in
%               valid calls of 'all_clean' and 'nzw_clean', respectively.
%
%     Phi     - the sparse tensor object from LiFE model to find voxel /
%               streamline node intersection.
%
%     dim     - a vector array of voxel size in mm, i.e. [ 1.25 1.25 1.25 ]
%
%     clobber - overwrite fields if they've already been computed
%
% OUTPUTS:
%     pconn - the paired connection object with the voxel indices for each 
%             edge added
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
% vx            = [ 1.25 1.25 1.25 ];
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % find unique voxels in each edge and add them to pconn
% pconn = fnFindPathVoxels(pconn, 'nzw', Phi, vx);
%
% Brent McPherson (c), 2017 - Indiana University
%

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

display('Finding voxels in tracts for all edges...')

% error without clobber set if volume is computed
if (isfield(pconn{1}.(label), 'volume') && clobber == 0)
    error('Volumes for this label are already computed. Set clobber = 1 to recompute.');
end    

% catch the total number of voxels and volume units
nvox = size(Phi, 2);
mms3 = prod(dim);

% estimate total wm vol
wmvol = nvox * mms3;

display([ 'Total white matter volume estimated at: ' num2str(wmvol) ' mm^3.' ])

tic;
parfor ii = 1:length(pconn)
    
    % pull the connection
    conn = pconn{ii}.(label);
    
    % if the connection is empty, fill in zero
    if isempty(conn.indices)
        
        % fill in empty voxel coords
        conn.volume.pvoxels = [];
    
    else
        
        % pull subtensor of the connection
        subtensor = Phi(:, :, conn.indices);
        
        % pull the unique voxel indices of the connection
        conn.volume.pvoxels = unique(subtensor.subs(:, 2));
    
    end
    
    % re-assign connection
    pconn{ii}.(label) = conn;
    
    % compute and store volume and proportional volume occupied by edge
    vol = size(conn.volume.pvoxels, 1) * mms3;
    pconn{ii}.(label).matrix.volume = vol;
    pconn{ii}.(label).matrix.prpvol = vol / wmvol;
        
end
time = toc;

display([ 'Found the voxels for all edges in ' num2str(time/60) ' minutes.' ]);

end

