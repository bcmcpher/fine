function [ pconn, tprof ] = feTractProfilePairedConnections(fe, pconn, label, msvol, mslab)
%feTractProfilePairedConnections runs virtual lesion on a field of pconn indices.
%
% INPUTS:
% - 'fe' is a fit fe structure that is used to build pconn
% - 'pconn' is the paired connection structure that stores fiber indices
%   between each unique pair of connections in the 
% - 'label' is the set of indices in pconn to perform virtual lesions for
% - 'msvol' is the aligned microstructural volume to extract tract profiles from
%
% OUTPUTS:
% - 'pconn' is the paired connections object with VL data stored and added
%   to the matrix field for generating networks
% - 'tprof' is the cell array of profiled output only as a cell array. For debugging.
%

% load msvol if it's not already loaded
if isstring(msvol)
    msvol = niftiRead(msvol);
end

display('Converting streamlines to ACPC space...');

% create fiber group
fg = feGet(fe, 'fg acpc');

% run parallelized tract profiles

% number of nodes to resample ms values along
nnodes = 100;

% check if average length is long enough for a profile to be reasonable?
% do I need a minimum number of streamlines to reasonable create a profile?

% preallocate tract profile output
tprof = cell(length(pconn), 1);
tpcnt = 0;

display(['Computing tract profiles on ' num2str(length(pconn)) ' connections...']);

tic;
for ii = 1:length(pconn)
    
    % pull field requested
    tmp = getfield(pconn{ii}, label);
    
    if size(tmp.indices, 1) > 2
        
        % create tract-wise fg
        tract = fgCreate('fibers', fg.fibers(tmp.indices));
        
        % compute tract profile
        tprof{ii} = dtiComputeDiffusionPropertiesAlongFG(tract, msvol, [], [], nnodes);
        
    else
        
        % skip empty connection
        tprof{ii} = [];
        continue
    end
    
end
time = toc;

display(['Computed ' num2str(tpcnt) ' tract profiles in ' num2str(round(time)/60) ' minutes.']);

clear ii time

% add virtual lesion to pconn object
for ii = 1:length(pconn)
        
    % pull subset field
    tmp = getfield(pconn{ii}, label);
    
    % look for an existing profile field
    if ~isempty(getfield(tmp, 'profile')) % if there is one
        
        % pull the existing fields and add the new one
        profile = getfield(tmp, 'profile');
        profile = setfield(mlab, tprof{ii});
        
    else
        
        % create the profile field and add the 
        profile = struct(mslab, tprof{ii});
    end
        
    % add whole vl output
    tmp = setfield(pconn{ii}, 'profile', profile);
       
    % reassign virtual lesions to paired connection array
    pconn{ii} = setfield(pconn{ii}, label, tmp);
    
end

end
