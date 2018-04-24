function [ pconn, tprof ] = feTractProfilePairedConnections(fg, pconn, label, msobj, mslab, nnodes, minNum, clobber)
%feTractProfilePairedConnections runs creates tract profile(s) for every existing 
% connection given a particular input. A nifti returns a single labeled
% profile, a dt6 structure will return 7 profiles in a new structure within
% the profile field.
%
% INPUTS:
%     fg      - fiber group in acpc space
%     pconn   - paired connection object created with fg
%     label   - string indicating the fiber groups for which to create tract profiles
%               either:
%                      'all' for all assigned streamlines or
%                      'nzw' for non-zero weighted fibers returned by LiFE
%     msobj   - a loaded dt6 or microstructural nifti image to compute profiles
%     mslab   - the label the new tract profile will be stored under
%               If a dt6 is passed, mslab is reset to 'dt6'
%     nnodes  - the number of nodes to which tract profiles will be sampled
%     minNum  - the minimum number of streamlines for a profile to be computed
%     clobber - overwrite existing profile label if it already exists
%             
% OUTPUTS:
%     pconn - is the paired connections object with the tract profile(s) added
%             in a field called 'profile'.
%
%     tprof - debugging output; the cell array that is added internally to pconn
%
% TODO:
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
% dt        = dtiloadDt6('dt6.mat');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create a profile for every connection with non-zero weighted fibers
% pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', favol, 'fa');
%
% % create a profile for every connection with non-zero weighted fibers
% pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', dt, 'dt6');
%
% Brent McPherson (c), 2017 - Indiana University
%

% parse optional arguments
if(~exist('nnodes', 'var') || isempty(nnodes))
    nnodes = 100;
end

if(~exist('minNum', 'var') || isempty(minNum))
    minNum = 3;
end

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

if(clobber == 1)
    display(['Overwriting ' mslab ' profiles...']);
end

if(isfield(msobj, 'dt6'))
    display('Computing tract profiles from dt6.');
    isdt6 = 1;
    mslab = 'dt6';
else
    display('Computing tract profiles from volume.');
    isdt6 = 0;
end

% check if the first profile label already exists alongside clobber for a
% faster exit from function if things shouldn't be recomputed.
if isfield(pconn{1}.(label), 'profile')
    if isfield(pconn{1}.(label).profile, mslab) && clobber == 0;
        error('Tract profiles with ''%s'' label already exist.\nPlease set clobber to 1 to explicitly overwrite existing profiles.', mslab);
    end
end

% check if average length is long enough for a profile to be reasonable?

% preallocate tract profile output
tprof = cell(length(pconn), 1);
tpcnt = 0;
tptry = 0;

disp(['Computing tract profiles on ' num2str(length(pconn)) ' connections...']);

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
            
            if(isdt6)
                
                [ tprof{ii}.dt6_fa, tprof{ii}.dt6_md, tprof{ii}.dt6_rd, ...
                  tprof{ii}.dt6_ad, tprof{ii}.dt6_cl, ~, ~, ...
                  tprof{ii}.dt6_cp, tprof{ii}.dt6_cs ] = ...
                  dtiComputeDiffusionPropertiesAlongFG(tract, msobj, [], [], nnodes);
                
            else
                
                % compute tract profile
                tprof{ii} = dtiComputeDiffusionPropertiesAlongFG(tract, msobj, [], [], nnodes);
                
            end
            
            % track how many connections are profiled
            tpcnt = tpcnt + 1;
            tptry = tptry + 1;
            
        catch
            
            warning(['Connection: ' num2str(ii) ' failed to compute profile.']);
            tptry = tptry + 1;
            
            if(isdt6)
                tprof{ii}.dt6_fa = nan(nnodes, 1);
                tprof{ii}.dt6_md = nan(nnodes, 1);
                tprof{ii}.dt6_rd = nan(nnodes, 1);
                tprof{ii}.dt6_ad = nan(nnodes, 1);
                tprof{ii}.dt6_cl = nan(nnodes, 1);
                tprof{ii}.dt6_cp = nan(nnodes, 1);
                tprof{ii}.dt6_cs = nan(nnodes, 1);
            else
                
                tprof{ii} = nan(nnodes, 1);
                
            end
        end
        
    else
        
        if(isdt6)
            tprof{ii}.dt6_fa = nan(nnodes, 1);
            tprof{ii}.dt6_md = nan(nnodes, 1);
            tprof{ii}.dt6_rd = nan(nnodes, 1);
            tprof{ii}.dt6_ad = nan(nnodes, 1);
            tprof{ii}.dt6_cl = nan(nnodes, 1);
            tprof{ii}.dt6_cp = nan(nnodes, 1);
            tprof{ii}.dt6_cs = nan(nnodes, 1);
            
        else
            
            % skip empty connection
            tprof{ii} = nan(nnodes, 1);
            continue
            
        end
        
    end
        
end
time = toc;

if(isdt6)
   tpcnt = tpcnt;
   tptry = tptry;
end

disp(['Computed ' num2str(tpcnt)  ' of ' num2str(tptry) ' possible tract profiles in ' num2str(round(time)/60) ' minutes.']);

clear ii time

%% add tract profile to pconn object

disp('Adding tract profiles to pconn...');

for ii = 1:length(pconn)
        
    % pull subset field
    tmp = pconn{ii}.(label);
    
    % look for an existing profile field
    if isfield(tmp, 'profile') % if there is one
        
        % pull the existing profiles and add the new one
        prof = tmp.profile;
        
        if(isdt6)
            prof.(mslab).fa = tprof{ii}.dt6_fa;
            prof.(mslab).md = tprof{ii}.dt6_md;
            prof.(mslab).rd = tprof{ii}.dt6_rd;
            prof.(mslab).ad = tprof{ii}.dt6_ad;
            prof.(mslab).cl = tprof{ii}.dt6_cl;
            prof.(mslab).cp = tprof{ii}.dt6_cp;
            prof.(mslab).cs = tprof{ii}.dt6_cs;
        else
            % add extra isfield check?
            prof.(mslab) = tprof{ii};
        end
        
    else
        
        if(isdt6)
            prof.(mslab).fa = tprof{ii}.dt6_fa;
            prof.(mslab).md = tprof{ii}.dt6_md;
            prof.(mslab).rd = tprof{ii}.dt6_rd;
            prof.(mslab).ad = tprof{ii}.dt6_ad;
            prof.(mslab).cl = tprof{ii}.dt6_cl;
            prof.(mslab).cp = tprof{ii}.dt6_cp;
            prof.(mslab).cs = tprof{ii}.dt6_cs;
        else
            
            % create the profile field and add the data
            prof = struct(mslab, tprof{ii});
        end
        
    end
        
    % add tract profile(s) to tmp
    tmp.profile = prof;
       
    % reassign tmp to a paired connection cell array
    pconn{ii}.(label) = tmp;
    
end

end
