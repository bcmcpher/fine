function [ pconn, tprof ] = feTractProfilePairedConnections(fg, pconn, label, msvol, mslab)
%feTractProfilePairedConnections runs virtual lesion on a field of pconn indices.
%
% INPUTS:
% - 'fg' is the optimized connectome returned by LiFE
%          fg = feGet(fe,'fibers acpc'); 
%           w = feGet(fe,'fiber weights');
%          fg = fgExtract(fg, 'keep', w > 0);
%
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

% load msvol if msvol is not already loaded
if isstring(msvol)
    msvol = niftiRead(msvol);
end

display('Converting streamlines to ACPC space...');

% run parallelized tract profiles

% number of nodes to resample ms values along
nnodes = 100;

% minimum number of streamlines
minNum = 4;

% check if average length is long enough for a profile to be reasonable?

% preallocate tract profile output
tprof = cell(length(pconn), 1);
tpcnt = 0;

display(['Computing tract profiles on ' num2str(length(pconn)) ' connections...']);

tic;    
fibers = fg.fibers;
parfor ii = 1:length(pconn)

    % pull field requested
    tmp = pconn{ii}.(label);
    
    % in testing, need minimum of 4 streamlines to compute profile
    if size(tmp.indices, 1) > minNum
        
        % create tract-wise fg
        
        tract = fgCreate('fibers', fibers(tmp.indices));
        
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

display(['Computed ' num2str(tpcnt) ' tract profiles in ' num2str(round(time)/60) ' minutes.']);

clear ii time

% %%%%
% %%%%  TODO: add logic to prevent overwriting an existing field %%%
% %%%%

% add virtual lesion to pconn object
parfor ii = 1:length(pconn)
        
    % pull subset field
    tmp = pconn{ii}.(label);
    
    % look for an existing profile field
    if isfield(tmp, 'profile') % if there is one
        
        % pull the existing profiles and add the new one
        prof          = tmp.profile;
        
        % add extra isfield check?
        prof.(mslab) = tprof{ii};
        
    else
        % create the profile field and add the data
        prof = struct(mslab, tprof{ii});
    end
        
    % add tract profile(s) to tmp
    tmp.profile = prof;
       
    % reassign tmp to a paired connection cell array
    pconn{ii}.(label) = tmp;
    
end

end
