function [ fh, cmap ] = fnPlotModuleAverages(tpMn, tpSd, mtMn, tpgrp, mtgrp, clr)
%fnPlotModuleAverages() creates the average profile and module matrix from
% a cell array of inputs. Recreates Figure ##.# from paper TBD.
%
%   Still in development
%

% TODO:
%
% high priority
% - generalize to N mtgrps passed - currenltly only handles up to 5
% - compute tract profile std if it's not passed
% - if 1 set of groupings are passed, apply to both profiles and matrices
% - - same argument based on dimensions?
%
% low priority
% - figure out lower / upper bounds of profile
% - be able to manually define the color scale?
% - multiple color scales?
% - multiple measures per group?

%% parse inputs for dimensions

% CHECK THAT ALL INPUT DIMENSIONS MATCH

% infer the number of subjects
%nsubj = size(grp, 1);

% pull the number of communities from labels
ncomm = 9;

% profile length / nodes
x = 1:100;

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
    
    % split data by groups
    gtpdat = tpMn(tpgrp == tpgrps(igrp));
    gtpstd = tpSd(tpgrp == tpgrps(igrp));
    
    % merge data into structes that can be summarized across
    tpmdat = cat(4, gtpdat{:});
    tpsdat = cat(4, gtpstd{:});
    
    % compute group-wise mean of data
    gMnTp{igrp} = nanmean(tpmdat, 4);
    gSeTp{igrp} = nanmean(tpsdat, 4);
        
end

clear igrp gtpdat gmndat tpdat mtdat

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

% % put fixed 2 groups into the data here
% 
% % feed created cell array into the specific outs here
% % then generalize the number of groups to be plotted
% hcMnTp = gMnTp{1};
% clMnTp = gMnTp{2};
% 
% hcSeTp = gSeTp{1};
% clSeTp = gSeTp{2};
% 
% hcMnMt = gMnMt{1};
% clMtMt = gMnMt{2};

% the number of groups needs to determine the patch sizes

% colors = {'red', 'blue', 'red', 'blue', 'red', 'blue'};

switch nmtgrp
    
    case 1
        % 1 group
        ptc = {[ 0 1 1 0 ], [ 1 1 0 0 ]};
        
    case 2
        % 2 group
        ptc = {[ 0 1 1 ], [ 1 1 0 ];  % upper
            [ 0 0 1 ], [ 1 0 0 ]}; % lower
        
    case 3
        % 3 group
        ptc = {[ 0.33 1 1 ], [ 1 1 0.33 ]; % upper
            [ 0 0 .33 1 1 .67 ], [ .67 1 1 .33 0 0]; % center
            [ 0 0 0.67], [ 0.67 0 0 ]}; % lower
        
    case 4
        % 4 group
        ptc = {[ 0.50 1 1 ],   [ 1 1 0.50 ];   % upper tri
            [ 0 0.50 1 1 ], [ 1 1 0.50 0 ]; % upper quad
            [ 0 0 1 0.50 ], [ 0.50 1 0 0 ]; % lower quad
            [ 0 0 0.50 ],   [ 0.50 0 0 ]};  % lower tri
        
    case 5
        % 5 group
        ptc = {[ 0.75 1 1 ], [ 1 1 0.75 ];         % upper tri
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
    %axis equal; axis square; axis tight;
    set(gca, 'XTick', [], 'YTick', []);
    
    % determine how to best fit y-axes?
    % set(gca, 'XLim', [ 0 size(hcMnTp, 3) ], ...
    %     'YLim', [ min([ ub1, ub2 ]) max([ lb1, lb2 ]) + .1 ]);
    
    hold off;
    
end

% define range
%cmap = hot(ceil(max([ clMtMt(:); hcMnMt(:) ])));
cmap = hot(ceil(max(max(max(gMnMt{:})))));

for ii = 1:size(pairs, 1)

    % create subplot
    subplot(ncomm, ncomm+1, upi(ii)); hold on;
    
    for imt = 1:nmtgrp
        
        % pull value from data
        val = round(gMnMt{imt}(pairs(ii, 1), pairs(ii, 2)));
        
        % pull colormap associated with the label
        color = cmap(val, :);
        
        % fill patch with scaled heatmap color
        patch(ptc{imt, :}, color, 'LineWidth', 1);
        
    end
        
    % drop all the extra panel labels
    axis equal; axis square; axis tight; 
    set(gca, 'XTick', [], 'YTick', []);
    
    hold off;
    
end

end

