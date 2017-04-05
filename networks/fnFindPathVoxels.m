function [ pconn ] = fnFindPathVoxels(Phi, pconn, label)
%fnFindPathVoxels finds the tensor indices of voxels for every non-zero
% connection to create either microstructural summaries of networks or for
% generating link networks.
%
% Phi is the sparse tensor - needed for finding the voxel indices
% Phi = feGet(fe, 'Phi');
%
%   This is hopefully a temporary solution as this functionality is useful
%   but is - apparently - very computationally demanding. I want this to be
%   found by default, but until it's fast I'll keep it separate.
%
%   Part of why it is slow is I had to compute it twice. Perhaps on demand
%   is always better, even if I basically always do it.
%

display('Finding voxels in tract for all edges...')

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

