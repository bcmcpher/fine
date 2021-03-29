function [ netw ] = fnTractProfileEdges(netw, fg, msobj, mslab, nnodes, minNum, clobber)
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
% dt        = dtiLoadDt6('dt6.mat');
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
    disp(['Overwriting ' mslab ' profiles...']);
end

if(isfield(msobj, 'dt6'))
    disp('Computing all possible tract profiles from dt6.');
    disp('Provided ''mslab'' will be ignored. dt6_* will prepend all new labels.');
    isdt6 = 1;
    mslab = 'dt6_fa';
else
    disp([ 'Computing tract profiles from ''' mslab ''' volume.' ]);
    isdt6 = 0;
end

% check if the first profile label already exists alongside clobber for a
% faster exit from function if things shouldn't be recomputed.
if isfield(netw.edges{1}, 'profile')
    disp('Some profiles have already been computed.');
    if isfield(netw.edges{1}.profile, mslab) && clobber == 0
        error('Tract profiles with ''%s'' label already exist.\nSet clobber to 1 to explicitly overwrite the pre-existing profile(s).', mslab);
    end
end

% check if the length is too short for a profile to be reasonable?

% preallocate tract profile output
tpcnt = 0;
tptry = 0;

% pull a quick count of edges w/ streamlines greater than the minimum
pc = sum(cellfun(@(x) size(x.fibers.indices, 1) > minNum, netw.edges));
disp(['Computing tract profiles on ' num2str(pc) ' present connections...']);

% pull lists for improved overhead in parallel loop
fibers = fg.fibers;
pconn = netw.edges;

tic;    
for ii = 1:size(pconn, 1)

    % pull next edge and roi indices
    edge = pconn{ii};
    r1_idx = netw.parc.pairs(ii, 1);
    r2_idx = netw.parc.pairs(ii, 2);
    
    % in testing, need minimum of 4 streamlines to compute profile
    if size(edge.fibers.indices, 1) > minNum
        
        % create tract-wise fg and reorient / resample
        tract = fgCreate('fibers', fibers(edge.fibers.indices));
        [ tfg, epi ] = dtiReorientFibers(tract, nnodes); % in theory don't have to resample
        
        % pull roi centers in acpc space
        roi1 = netw.nodes{r1_idx}.center.acpc';
        roi2 = netw.nodes{r2_idx}.center.acpc';
        
        % find the distance between the start of the profile and each roi center
        epi_roi1 = norm(epi - roi1);
        epi_roi2 = norm(epi - roi2);
        
        % if roi2 is closer to the start of the tract profile than roi1, lrflip the fibers
        if (epi_roi2 < epi_roi1)
            tfg.fibers = cellfun(@(x) fliplr(x), tfg.fibers, 'UniformOutput', false);
        end
        
        % create superfiber representation for the edge on the flipped fibers
        if ~isfield(edge, 'superfiber')
            sf_name = strcat(netw.nodes{r1_idx}.name, '-', netw.nodes{r2_idx}.name);
            edge.superfiber = dtiComputeSuperFiberRepresentation(tfg, [], nnodes);
            edge.superfiber.name = sf_name;
        end
        
        % if the variance of the streamlines distance is far, too many
        % streamlines are dropped to reliably compute profile, so try / catch
        try
            
            if(isdt6)
                
                % compute all the tract profiles for a dt6
                [ edge.profile.dt6_fa, edge.profile.dt6_md, edge.profile.dt6_rd, ...
                  edge.profile.dt6_ad, edge.profile.dt6_cl, ~, ~, ...
                  edge.profile.dt6_cp, edge.profile.dt6_cs ] = ...
                  dtiComputeDiffusionPropertiesAlongFG(tfg, msobj, [], [], nnodes);

            else
                
                % compute the tract profile for the passed volume
                [ edge.profile.(mslab), ~, ~, ~, ~, ...
                  ~, ~, ~, ~ ]= dtiComputeDiffusionPropertiesAlongFG(tract, msobj, [], [], nnodes);

            end
            
            % track how many connections are profiled
            tpcnt = tpcnt + 1;
            tptry = tptry + 1;
            
        catch
            
            warning(['Connection: ' num2str(ii) ' failed to compute profile.']);
            tptry = tptry + 1;
            
            if(isdt6)
                edge.profile.dt6_fa = nan(nnodes, 1);
                edge.profile.dt6_md = nan(nnodes, 1);
                edge.profile.dt6_rd = nan(nnodes, 1);
                edge.profile.dt6_ad = nan(nnodes, 1);
                edge.profile.dt6_cl = nan(nnodes, 1);
                edge.profile.dt6_cp = nan(nnodes, 1);
                edge.profile.dt6_cs = nan(nnodes, 1);
            else
                edge.profile.(mslab) = nan(nnodes, 1);                
            end
        end
        
    else
        
        % too few streamlines exist to compute a profile / it's empty
        if ~isempty(edge.fibers.indices)
            warning([ 'Edge ' num2str(ii) ' is not empty: ' num2str(length(edge.fibers.indices)) ' streamline; less than ' num2str(minNum) ]);
            % zero this connection if it's this small?
            % edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
            % there's a bunch of other potentially non-zero fields though...
        end
        
        % fill in empty dt6 fields
        if(isdt6)
            edge.profile.dt6_fa = [];
            edge.profile.dt6_md = [];
            edge.profile.dt6_rd = [];
            edge.profile.dt6_ad = [];
            edge.profile.dt6_cl = [];
            edge.profile.dt6_cp = [];
            edge.profile.dt6_cs = [];            
        else
            % fill in empty single node profile
            edge.profile.(mslab) = [];
        end
        
    end
    
    % reassign the edge to pconn once profiles are computed
    pconn{ii} = edge;
        
end
time = toc;

% reassign connections
netw.edges = pconn;

disp(['Computed ' num2str(tpcnt)  ' of ' num2str(tptry) ' possible tract profiles in ' num2str(round(time)/60) ' minutes.']);

end
