function [ fh ] = fnPlotModuleGroups(mats, srt, clr)
%fnPlotModuleGroups() plots each unique module of each group as a point
% Points are sorted by within / between module connections in the order of
% the upper diagonal. Error bars are the standard error of the group.
%
% Recreates Figure ##.# from paper TBD.
%
% INPUTS:
%     mats - a N x N x subj matrix of network modules
%     srt  - a subj x 1 array grouping subjects together into categories
%     clr  - a ngrp x 1 cell array of colors to plot subjects' points
%
% OUTPUTS:
%     fh   - figure handle of the resulting plot
%
% Brent McPherson, (c) 2018 Indiana University
%

% find the unique number of group labels
grps = unique(srt);
ngrp = size(grps, 1);

% if colors / groups don't match, quit
if ngrp ~= size(clr, 1)
    error('The number of groups does not match the number of color labels.');
end

% pull the number of modules from input
nmod = size(mats, 1);

% unique group should match the # of colors passed
nsubj = size(mats, 3);

% preallocate output data
gmn = nan(nmod, nmod, ngrp);
gse = nan(nmod, nmod, ngrp);

% compute average / se for each group
for grp = 1:ngrp
    
    % pull label
    grpl = grps(grp);
    
    % pull nsubj in grp
    tsbj = sum(srt == grpl);
    
    % pull data
    tmat = mats(:, :, srt == grpl);
    
    % compute the mean / se of the group
    gmn(:,:,grp) = mean(tmat, 3, 'omitnan');
    gse(:,:,grp) = std(tmat, [], 3, 'omitnan') ./ sqrt(tsbj);

end

clear grp grpl tsbj tmat

% pull maximum range of axes
amax = max(mats(:));
amid = amax / 2;

% create diagonal indices from input size
diag = [ 1:nmod; 1:nmod ]';

% create the upper diagonal combinations
updg = nchoosek(1:nmod, 2);

% create the index order to plot all the points
indx = [ diag; updg ];

% total number of modules (data points) to plot per group
npts = size(indx, 1);

% create module labels for axis
labs = cell(npts, 1);
for lab = 1:npts
    labs{lab} = [ num2str(indx(lab, 1)) '-' num2str(indx(lab, 2)) ];
end

% define jitter between points
jitter = linspace(-0.20, 0.20, length(1:ngrp));

% plot the data
fh = figure('Position', [ 150 575 1675 625 ]); hold on;

% for every subject
for grp = 1:ngrp
    
    % pull subject data and color from input
    mnmat = gmn(:, :, grp);
    semat = gse(:, :, grp);
    color = clr{grp};
    
    % for every module
    for mod = 1:size(indx, 1)
        
        % the data value to plot
        pt = mnmat(indx(mod, 1), indx(mod, 2));
        eb = semat(indx(mod, 1), indx(mod, 2));
      
        % plot the data point
        plot(mod + jitter(grp), pt, 'o', 'MarkerFaceColor', color, ...
            'MarkerEdgeColor', 'black', 'LineWidth', 0.25, 'MarkerSize', 8);
        
        % plot the error bar
        plot([ mod + jitter(grp) mod + jitter(grp) ], [ pt - eb pt + eb ], ...
             'black', 'LineWidth', 0.25);
        
    end
    
end

% format the plot axes
set(gca, 'YTick', [ 0 amid amax ], 'YTickLabel', [ 0 amid amax ], 'YLim', [ 0 amax + (amax * .05) ]);
set(gca, 'XTick', 1:npts, 'XTickLabel', labs, 'XLim', [ 0 npts + 1 ], 'TickLength', [ 0 0 ]);
set(gca, 'XTickLabelRotation', 90);

end

