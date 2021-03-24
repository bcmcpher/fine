function [ netw, wmvol ] = fnFindPathVoxels(netw, Phi, clobber)
%%
% REPLACED WITH fnAverageEdgePropertyS() TO NOT REQUIRE ENCODE STRUCTURE
% AND FOLLOW UP WITH fnAverageEdgeProperty()
%%
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
%     wmvol - the total white matter volume according to the provided voxel 'dim'
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

disp('Finding voxels in tracts for all edges...')

% error without clobber set if volume is computed
if (isfield(netw.pconn{1}, 'volume') && clobber == 0)
    error('Volumes for this label are already computed. Set clobber = 1 to recompute.');
end    

% catch the total number of voxels and volume units
netw.parc.wmvox = size(Phi, 2);
mms3 = prod(netw.parc.voxmm);

% estimate total wm vol
wmvol = netw.parc.wmvox * mms3;

disp([ 'Total white matter volume estimated at: ' num2str(wmvol) ' mm^3.' ])

% store computed white matter volume
netw.parc.wmvol = wmvol;

tic;
for ii = 1:length(netw.pconn)
    
    % pull the connection
    conn = netw.pconn{ii};
    
    % if the connection is empty, fill in zero
    if isempty(conn.fibers.indices)
        
        % fill in empty voxel coords
        conn.pvoxels = [];
    
    else
        
        % pull subtensor of the connection
        subtensor = Phi(:, :, conn.fibers.indices);
        
        % pull the unique voxel indices of the connection
        conn.pvoxels = unique(subtensor.subs(:, 2));
    
    end
        
    % compute and store volume and proportional volume occupied by edge
    vol = size(conn.pvoxels, 1) * mms3;
    conn.volume.volume = vol;
    conn.volume.prpvol = vol / wmvol;
    
    % re-assign connection
    netw.pconn{ii} = conn;
    
end
time = toc;

disp([ 'Found the voxels for all edges in ' num2str(time/60) ' minutes.' ]);

end

