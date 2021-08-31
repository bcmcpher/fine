function [ netw ] = fnTractProfileEdges(netw, fg, msobj, mslab, nnodes, minNum, clobber)
%feTractProfilePairedConnections runs creates tract profile(s) for every existing 
% connection given a particular input. A nifti returns a single labeled
% profile, a dt6 structure will return 7 profiles in a new structure within
% the profile field.
%
% INPUTS:
%     netw    - fiber group in acpc space
%     fg      - paired connection object created with fg
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
    clobber = false;
end

if(clobber == 1)
    disp('''clobber'' set to ''true''. Potentially overwriting stored profiles...');
end

% check if the first profile label already exists alongside clobber for a
% faster exit from function if things shouldn't be recomputed.
if isfield(netw.edges{1}, 'profile')
    disp('Some profiles have already been computed.');
    prevProf = intersect(fieldnames(netw.edges{1}.profile), mslab);
    if ~isempty(prevProf) && ~clobber
        fprintf('Tract profiles with the following labels already exist:\n');
        fprintf(' - %s\n', prevProf{:});
        error('Set ''clobber'' to true to overwrite currently stored profiles.');
    end
end

% check if there are cell array of input names / volumes
if iscell(mslab) || iscell(msobj)
    % check if they're the same length
    if (length(mslab) == length(msobj))
        % if they labels / volumes have the same length, check dimensions
        dims = cellfun(@(x) size(x.data), msobj, 'UniformOutput', false);
        if isequal(dims{:})
            disp('Array of inputs have matching dimensions.');
        else
            warning('Dimensions of input volumes do not match. This may cause issues.');
        end
        ninputs = length(mslab);
    else
        error('The cell array of labels and volumes are not the same length. Unable to appropriately assign inputs.');
    end
else
    
    % make the single inputs cells for indexing sanity
    mslab = {mslab};
    msobj = {msobj};
    
end

if(isfield(msobj, 'dt6'))
    disp('Computing all possible tract profiles from dt6.');
    disp('dt6_* will prepend all new labels.');
    isdt6 = 1;
    mslab = {'dt6_fa'};
else
    %disp([ 'Computing tract profiles from ''' mslab ''' volume.' ]);
    fprintf('Computing tract profiles for all edges of:\n');
    fprintf(' - %s\n', mslab{:});

    isdt6 = 0;
end

% check if the length is too short for a profile to be reasonable?

% preallocate tract profile output
tpcnt = 0;
tptry = 0;

% pull a quick count of edges w/ streamlines greater than the minimum
pc = sum(cellfun(@(x) size(x.fibers.indices, 1) > minNum, netw.edges));
disp(['Computing tract profiles on ' num2str(pc) ' present connections...']);

% pull lists for improved overhead in loop
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
        % in theory don't have to resample
        [ tfg, epi ] = dtiReorientFibers(tract, nnodes); % in practice you do
        
        % pull roi centers in acpc space
        roi1 = netw.nodes{r1_idx}.center.acpc';
        roi2 = netw.nodes{r2_idx}.center.acpc';
        
        % find the distance between the start of the resampled profile and each roi center
        epi_roi1 = norm(epi - roi1);
        epi_roi2 = norm(epi - roi2);
        
        % if roi2 is closer to the start of the tract profile than roi1, lrflip the fibers
        if (epi_roi2 < epi_roi1)
            tfg.fibers = cellfun(@(x) fliplr(x), tfg.fibers, 'UniformOutput', false);
        end
        % profiles should all be oriented / stored in i -> j node orientation (upper diagonal)
        
        % create superfiber representation for the edge on the flipped fibers
        if ~isfield(edge, 'superfiber')
            sf_name = strcat(netw.nodes{r1_idx}.name, '-', netw.nodes{r2_idx}.name);
            edge.superfiber = dtiComputeSuperFiberRepresentation(tfg, [], nnodes);
            edge.superfiber.name = sf_name;
        end
        
        % create the center of the average profile
        cfib = cat(3, tfg.fibers{:});
        edge.profile.center = mean(cfib, 3)';
        
        % if the variance of the streamlines distance is far, too many
        % streamlines are dropped to reliably compute profile, so try / catch
        
        % for every volume
        %for jj = 1:ninputs
        
        % type to compute the profile on the oriented reampled edge
        try
            
            for jj = 1:ninputs
                
                % if it's a dt6
                if(isdt6)
                    
                    % compute all the tract profiles for a dt6
                    [ edge.profile.dt6_fa, edge.profile.dt6_md, edge.profile.dt6_rd, ...
                        edge.profile.dt6_ad, edge.profile.dt6_cl, ~, ~, ...
                        edge.profile.dt6_cp, edge.profile.dt6_cs ] = ...
                        dtiComputeDiffusionPropertiesAlongFG(tfg, msobj{jj}, [], [], nnodes);
                    
                else
                    
                    % compute the tract profile for the passed volume
                    [ edge.profile.(mslab{jj}), ~, ~, ~, ~, ...
                        ~, ~, ~, ~ ]= dtiComputeDiffusionPropertiesAlongFG(tract, msobj{jj}, [], [], nnodes);
                    
                end
            end
            
            % track how many connections are profiled
            tpcnt = tpcnt + 1;
            tptry = tptry + 1;
            
        catch
            
            % if 1 fails all are zeroed (?)
            warning(['Connection: ' num2str(ii) ' failed to compute profile.']);
            tptry = tptry + 1;
            
            for jj = 1:ninputs
                if(isdt6)
                    edge.profile.dt6_fa = nan(nnodes, 1);
                    edge.profile.dt6_md = nan(nnodes, 1);
                    edge.profile.dt6_rd = nan(nnodes, 1);
                    edge.profile.dt6_ad = nan(nnodes, 1);
                    edge.profile.dt6_cl = nan(nnodes, 1);
                    edge.profile.dt6_cp = nan(nnodes, 1);
                    edge.profile.dt6_cs = nan(nnodes, 1);
                else
                    edge.profile.(mslab{jj}) = nan(nnodes, 1);
                end
            end
        end
        %end
        
    else
        
        % too few streamlines exist to compute a profile / it's empty
        if ~isempty(edge.fibers.indices)
            warning([ 'Edge ' num2str(ii) ' is not empty: ' num2str(length(edge.fibers.indices)) ' streamline; less than ' num2str(minNum) ]);
            % zero this connection if it's this small?
            % edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
            % there's a bunch of other potentially non-zero fields though...
        end
        
        % fill in an empty average profile, even if it's not truly empty
        edge.profile.center = [];
        
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
            
            % fill in every value empty (?) - If one fails all get zeroed?
            for jj = 1:ninputs
                % fill in empty single node profile
                edge.profile.(mslab{jj}) = [];
            end
            
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
