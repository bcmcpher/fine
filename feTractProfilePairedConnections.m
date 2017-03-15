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

% load msvol if it'smsv not already loaded
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
parfor ii = 1:length(pconn)
    
    % pull field requested
    tmp = getfield(pconn{ii}, label);
    
    % in testing, need minimum of 4 streamlines to compute profile
    if size(tmp.indices, 1) > 3
        
        % create tract-wise fg
        tract = fgCreate('fibers', fg.fibers(tmp.indices));
        
        % compute tract profile
        tprof{ii} = dtiComputeDiffusionPropertiesAlongFG(tract, msvol, [], [], nnodes);
        
        % track how many connections are profiled
        tpcnt = tpcnt + 1;
        
    else
        
        % skip empty connection
        tprof{ii} = [];
        continue
    end
        
end
time = toc;

display(['Computed ' num2str(tpcnt) ' tract profiles in ' num2str(round(time)) ' seconds.']);

clear ii time

% TODO: add logic to prevent overwriting an existing field

% add virtual lesion to pconn object
for ii = 1:length(pconn)
        
    % pull subset field
    tmp = getfield(pconn{ii}, label);
    
    % look for an existing profile field
    % NEED TO FIX SO WHEN IT'S FOUND, THIS IS TRUE
    if isempty(strfind(fieldnames(tmp), 'profile')) % if there is one
        
        % pull the existing fields and add the new one
        prof = getfield(tmp, 'profile');
        prof = setfield(mlab, tprof{ii});
        
    else
        
        % create the profile field and add the data
        prof = struct(mslab, tprof{ii});
    end
        
    % add tract profile(s) to tmp
    tmp = setfield(tmp, 'profile', prof);
       
    % reassign tmp to a paired connection cell array
    pconn{ii} = setfield(pconn{ii}, label, tmp);
    
end

end