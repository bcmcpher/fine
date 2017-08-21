function [ pconn, tprof ] = feTractProfilePairedConnections(fg, pconn, label, msvol, mslab)
%feTractProfilePairedConnections runs creates tract profile(s) for every existing 
% connection given a particular input.
%
% INPUTS:
%     fg    - fiber group in acpc space
%     pconn - paired connection object created with fg
%     label - string indicating the fiber groups for which to create tract profiles
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%     msvol - a loaded microstructural nifti image to compute profiles with
%     mslab - the label the new tract profile will be stored under
%             
% OUTPUTS:
%     pconn - is the paired connections object with the tract profile(s) added
%             in a field called 'profile'.
%
%     tprof - debugging output; the cell array that is added internally to pconn
%
% TODO:
% - pass dt6 instead of msvol to create multiple profiles
% - check if a labels exists and don't overwrite / add clobber
% - pass nnodes / minnum as fxn level arguments
% - check if the connections are too short for a reasonable profile (?)
%
% EXAMPLE:
%
% % load data
% parc  = niftiRead('labels.nii.gz');
% favol = niftiRead('fa.nii.gz');
% fg        = feGet(fe, 'fibers acpc');
% fibers    = fg.fibers;
% fibLength = fefgGet(fg, 'length');
% weights   = feGet(fe, 'fiberweights');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create a profile for every connection with non-zero weighted fibers
% pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', favol, 'fa');
%
% Brent McPherson (c), 2017 - Indiana University
%

% load msvol if msvol is not already loaded
if isstring(msvol)
    msvol = niftiRead(msvol);
end

% run parallelized tract profiles

% number of nodes to resample ms values along
nnodes = 100;

% minimum number of streamlines
minNum = 4;

% check if average length is long enough for a profile to be reasonable?

% preallocate tract profile output
tprof = cell(length(pconn), 1);
tpcnt = 0;
tptry = 0;

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
        
        % if the variance of the streamlines distance is far, too many
        % streamlines are dropped to reliably compute profile
        try
            
            % compute tract profile
            tprof{ii} = dtiComputeDiffusionPropertiesAlongFG(tract, msvol, [], [], nnodes);
            
            % track how many connections are profiled
            tpcnt = tpcnt + 1;
            tptry = tptry + 1;
        catch
            
            warning(['Connection: ' num2str(ii) ' failed to compute profile.']);
            tprof{ii} = nan(nnodes, 1);
            tptry = tptry + 1;
        end
        
    else
        
        % skip empty connection
        tprof{ii} = [];
        continue
    end
        
end
time = toc;

display(['Computed ' num2str(tpcnt)  ' of ' num2str(tptry) ' possible tract profiles in ' num2str(round(time)/60) ' minutes.']);

clear ii time

%% add virtual lesion to pconn object

parfor ii = 1:length(pconn)
        
    % pull subset field
    tmp = pconn{ii}.(label);
    
    % look for an existing profile field
    if isfield(tmp, 'profile') % if there is one
        
        % pull the existing profiles and add the new one
        prof = tmp.profile;
        
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
