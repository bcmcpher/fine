function [ fh, cmap ] = fnPlotModuleAverages(tpMn, tpSd, mtMn, tpgrp, mtgrp, clr, ilim)
%fnPlotModuleAverages() creates the average profile and module matrix from
% a cell array of inputs. The number of profiles and matrix modules may be 
% different. Depending on the data that is combined, this plot can display
% multiple modality observations within a subject or different groups of
% multiple subjects.
%
% Recreates Figure ##.# from paper TBD.
%
% INPUTS:
%     tpMn  - a 1xN cell array of mean tract profile tensors
%     tpSd  - a 1xN cell array of std tract profile tensors 
%             (optional) if this is passed empty, std is computed between tpMn observations
%     mtMn  - a 1xM cell array of network modules connectivity
%     tpgrp - a 1xN array assigning group membership for averaging profile inputs
%     mtgrp - a 1xM array assigning group membership for averaging module inputs
%     clr   - colors to plot group-wise profiles
%     ilim  - user specified y-axis limit range on average profile plots (optional)
%
% OUTPUTS:
%     fh    - figure handle of the resulting plot
%     cmap  - color map produced internally to map average network modules
%
% TODO:
%
% high priority
% - make sure color map indices increase correctly
% - generalize to N mtgrps passed - currenltly only handles up to 5
%
% low priority
% - be able to manually define the color scale?
% - multiple color scales?
%
% - multiple measures per group? - It already does, stack repeated measures
%   and group accordingly (grp1-fa, grp2-fa, grp1-md, grp2-md, etc)
%
% Brent McPherson, (c) 2018 Indiana University
%

%% parse inputs for dimension checks

% grab the dimensions of the input profiles
[ x, y, z ] = cellfun(@size, tpMn, 'UniformOutput', true);

% if any row / col is different or if there are different lengths
if ~all([ all(x) all(y) ]) || ~all(z)
    error('An entry in the individual profile means does not have the correct size. Dimensions must match.');
end

clear x y z

% check standard deviation if tpSd
if ~isempty(tpSd)
    
    % grab the dimensions of the input profiles
    [ x, y, z ] = cellfun(@size, tpSd, 'UniformOutput', true);
    
    % if any row / col is different or if there are different lengths
    if ~all([ all(x) all(y) ]) || ~all(z)
        error('An entry in the individual profile sd''s does not have the correct size. Dimensions must match.');
    end
    
    clear x y z
    
end

% for every input module matrix, get the row / col dim
[ x, y ] = cellfun(@size, mtMn, 'UniformOutput', true);

% if anything isn't the same dimension, error out
if ~all([ all(x) all(y) ])
    error('An entry in the individual modules does not have the correct size. Dimensions must match.');
end

clear x y

%% parse inputs for processing

% infer the number of subjects
nsubj = size(tpMn, 1);

% pull the number of communities from labels
ncomm = size(tpMn{1}, 1);

% profile length / nodes
nnode = size(tpMn{1}, 3);
x = 1:nnode;

if(~exist('ilim', 'var') || isempty(ilim))
    ilim = [];
    plim = 'auto';
end

% if something is passed
if ~isempty(ilim)
    
    % pass an array forward
    if isnumeric(ilim)
        plim = 'pass';
        limit = ilim;
    end
    
    % pass a string as what's requested
    if ischar(ilim)
        plim = ilim;        
    end
    
end

%% parse together input profile data based into defined group averages

% pull the number of groups and their count
tpgrps = unique(tpgrp);
ntpgrp = size(tpgrps, 1);
%ctpgrp = sum(bsxfun(@eq, tpgrp, tpgrps'))';

% preallocate group-wise cell array
gMnTp = cell(ntpgrp, 1);
gSeTp = cell(ntpgrp, 1);

% split inputs by group
for igrp = 1:ntpgrp
    
    % split data by groups, recombine to compute mean
    gtpdat = tpMn(tpgrp == tpgrps(igrp));
    tpmdat = cat(4, gtpdat{:});
    gMnTp{igrp} = mean(tpmdat, 4, 'omitnan');
    
    % if std is not passed, compute it; otherwise average what's passed
    if isempty(tpSd)
        
        % compute standard deviation
        gSeTp{igrp} = std(tpmdat, 0, 4, 'omitnan');
    
    else
        
        % split data by groups, recombine to compute mean
        gtpstd = tpSd(tpgrp == tpgrps(igrp));
        tpsdat = cat(4, gtpstd{:});
        gSeTp{igrp} = mean(tpsdat, 4, 'omitnan');
    
    end
        
end

clear igrp gtpdat tpmdat gmndat gtpstd tpsdat mtdat

%% extract min / max of vector for axes limits

switch plim
    
    case {'free'}
        % don't fix plot axes
        flim = false;
        
    case{'auto'}
        % find and apply min/max +- 5% limit
        
        % pull into array
        tpmdat = cat(4, tpMn{:});
        
        % better ub/lb computation? currently sucks...
        %
        % % grab global mean / std
        % tpmn = mean(tpmdat, 4, 'omitnan');
        % tpse = std(tpmdat, 0, 4, 'omitnan');
        %
        % % grap 3 std upper bound
        % tpub = tpmn + 3*tpse;
        %
        % pad = range(tpub) * 0.5;
        %
        % % create limit to pass forward
        % mnmx = minmax(tpub(:)');
        %
        % % create limit
        % limit = [ mnmx(1) - pad mnmx(2) + pad ];
        
        % reshape for finding min / max
        bndat = reshape(tpmdat, [ ncomm*ncomm*nsubj nnode ]);
        
        % grab the absolute upper / lower boundaries
        mnmx = minmax(bndat(:)');
        
        % pad min / max 5% of the total range
        pad = range(mnmx) * .05;
        limit = [ mnmx(1) - pad mnmx(2) + pad ];
        flim = true;
        
    case 'pass'
        
        flim = true;
        
    otherwise
        
        error('Unacceptable limits for profiles passed. Please pass a [ # # ] array.');
        
end

%% parse together input matrix data based into defined group averages

% pull the number of groups and their count
mtgrps = unique(mtgrp);
nmtgrp = size(mtgrps, 1);
%cmtgrp = sum(bsxfun(@eq, mtgrp, mtgrps'))';

% preallocate group-wise cell array
gMnMt = cell(nmtgrp, 1);

% split inputs by group
for igrp = 1:nmtgrp
    
    % split data by groups
    gmndat = mtMn(mtgrp == mtgrps(igrp));
    
    % merge data into structes that can be summarized across
    mtmdat = cat(3, gmndat{:});
    
    % compute group-wise mean of data
    gMnMt{igrp} = nanmean(mtmdat, 3);
        
end

clear igrp gmndat mtdat

% get the divisions of upper diagonal cells
switch nmtgrp
    
    case 1
        % 1 group
        ptc = {[ 0 1 1 0 ], [ 1 1 0 0 ]};
        
    case 2
        % 2 group
        ptc = {[ 0 1 1 ], [ 1 1 0 ];  % upper
            [ 0 0 1 ], [ 1 0 0 ]};    % lower
        
    case 3
        % 3 group
        ptc = {[ 0.33 1 1 ], [ 1 1 0.33 ]; % upper
            [ 0 0 .33 1 1 .67 ], [ .67 1 1 .33 0 0]; % center
            [ 0 0 0.67], [ 0.67 0 0 ]}; % lower
        
    case 4
        % 4 group
        ptc = {[ 0.50 1 1 ],   [ 1 1 0.50 ]; % upper tri
            [ 0 0.50 1 1 ], [ 1 1 0.50 0 ];  % upper quad
            [ 0 0 1 0.50 ], [ 0.50 1 0 0 ];  % lower quad
            [ 0 0 0.50 ],   [ 0.50 0 0 ]};   % lower tri
        
    case 5
        % 5 group
        ptc = {[ 0.75 1 1 ], [ 1 1 0.75 ];      % upper tri
            [ 0.25 .75 1 1 ], [ 1 1 0.75 .25];  % upper quad
            [ 0 0 .25 .25 1 1 .75 ], [ .75 1 1 1 .25 0 0 ]; % center
            [ 0 0 .75 0.25 ], [ 0.25 .75 0 0 ]; % lower quad
            [ 0 0 0.25 ], [ 0.25 0 0 ]};        % lower tri
        
    otherwise
        error('This function cannont currently support %s matrix groups.', nmtgrp);
end

% figure; hold on;
% patch(sqr{1, :}, 'green');
% set(gca, 'XLim', [ -.5 1.5 ], 'YLim', [ -.5 1.5]);
% axis square;
% 
% patch(ptc{1, :}, colors{1});
% patch(ptc{2, :}, colors{2});
% patch(ptc{3, :}, colors{3});
% patch(ptc{4, :}, colors{4});
% patch(ptc{5, :}, colors{5});
   
%% make the plot

% find the index of each panel of the figure
figi = reshape(1:(ncomm * (ncomm + 1)), [ ncomm + 1 ncomm ])';

% find the lower and upper diagonal indices
lpi = tril(figi);
lpi = lpi(lpi > 0);
upi = triu(figi, 1);
upi = sort(upi(upi > 0));

% find lower diagonal (between modules)
lowd = nchoosek(1:ncomm, 2);

% find diagonal (within modules)
diag = [ 1:ncomm; 1:ncomm ]';

% combine the module indices to compare
pairs = [ diag; lowd ];

% sort into upper diagonal indices
pairs = sortrows(pairs, [ 1 2 ]);

% loop and create all the module average (or overlaid) profiles
fh = figure('Position', [ 375 -225 1500 1125 ]); 

for ii = 1:size(pairs, 1)
    
    % create subplot
    subplot(ncomm, ncomm+1, lpi(ii)); hold on;
    
    % plot all the shaded regions
    for itp = 1:ntpgrp
        
        % grab individual profile
        mn = squeeze(gMnTp{itp}(pairs(ii, 1), pairs(ii, 2), :));
        se = squeeze(gSeTp{itp}(pairs(ii, 1), pairs(ii, 2), :));
        
        % grab upper and lower bounds
        ub = (mn - se)';
        lb = (mn + se)';
        
        % make the patch
        yP = [ lb, fliplr(ub) ];
        yP(isnan(yP)) = [];
        xP = [ x, fliplr(x) ];
        
        % fill in the shaded region between standard deviations
        patch(xP, yP, 1, 'facecolor', clr{itp}, 'edgecolor', 'none', 'facealpha', 0.20);
        
    end
    
    % plot all the lines on top of the shaded regions
    for itp = 1:ntpgrp
        
        % grab individual profile
        mn = squeeze(gMnTp{itp}(pairs(ii, 1), pairs(ii, 2), :));
        
        % plot the line
        plot(x, mn, 'color', clr{itp});
        
    end
    
    % drop all the extra panel labels
    % axis equal; axis square; axis tight;
    axis tight;
    set(gca, 'XTick', [], 'YTick', []);
    
    % determine how to best fit y-axes?
    set(gca, 'XLim', [ 0 nnode ])
    
    if flim
        set(gca, 'YLim', [ limit(1) limit(2) ]);
    end
    
    hold off;
    
end

%% define color mapping 

% find the min and max values in the matrix
crng = [ min(min(cat(2, gMnMt{:}))) max(max(cat(2, gMnMt{:}))) ];

% define the color map
cmap = hot(64);
caxis([ crng(1) crng(2) ]);

% create the spacing between min / max that corresponds to colors
vals = linspace(crng(1), crng(2), 64);

for ii = 1:size(pairs, 1)

    % create subplot
    subplot(ncomm, ncomm+1, upi(ii)); hold on;
    
    for imt = 1:nmtgrp
        
        % pull value from data
        val = gMnMt{imt}(pairs(ii, 1), pairs(ii, 2));
        
        % pull the index for the color map
        [ ~, cmv ] = min(abs(val - vals));
        
        % pull colormap associated with the label
        color = cmap(cmv, :);
        
        % fill patch with scaled heatmap color
        patch(ptc{imt, :}, color, 'LineWidth', 1);
        
    end
        
    % drop all the extra panel labels
    axis equal; axis square; axis tight; 
    set(gca, 'XTick', [], 'YTick', []);
    
    hold off;
    
end

end

