function [ fig ] = plotVirtualLesion(se, plt)
%plotVirtualLesion recreates the distributions from the LiFE paper for any
% virtual lesion output.
%
% INPUTS:
%     se  - an output from virtual lesion containing the data to plot.
%           Typically found in pconn{#}.(label).vl
%
%     plt - string indicating the specific plot to make
%             either:
%                     'all' (default) print all 3 plots
%                     '1' the RMSE probability distributions
%                     '2' the strength of evidence distributions
%                     '3' the distribution measures computed by virtual lesion
%
% OUTPUTS:
%     fig - a figure handle containing the axes of each requested figure
%
% TODO:
% - improve the defualt axes of the virtual lesion values
%
% EXAMPLE:
%
% % create plots from a specific edge
% fig = fnCreateLinkNetwork(pconn{2}.nzw.vl);
%
% % create a single plot from a specific edge
% fig = fnCreateLinkNetwork(pconn{2}.nzw.vl, 2);
%
% Brent McPherson (c), 2017 - Indiana University
%

% by default, plot all the figures
if(~exist('plt', 'var') || isempty(plt))
    plt = 'all';
else
    % force to a string
    plt = num2str(plt);
end

% figure index for they don't skip
fi = 1;

%% raw RMSE distirbutions

if (strcmp(plt, '1') || strcmp(plt, 'all'))
    
    % catch output
    fig(fi).name = 'rmse_distributions';
    fig(fi).h    = figure('name', fig(1).name, 'color', 'w');
    fig(fi).type = 'eps';
    
    % initialize header
    set(fig(fi).h, 'Units', 'normalized', 'Position', [ 0.005 0.50 0.30 0.35 ]);
    
    % better relative x-boundary
    dat = [ se.nolesion.rmse.all se.lesion.rmse.all ];
    jitter = (min(dat) * 0.05);
    mindat = round(min(dat) - jitter);
    maxdat = round(max(dat) + jitter);
    middat = round(mindat + ((maxdat - mindat) / 2));
    
    % better relative y-boundary
    ymax = max([ max(se.lesion.hist) max(se.nolesion.hist) ]);
    ymid = ymax / 2;
    yjit = ymax * 0.10;
    
    % plot the distributions
    plot(se.lesion.xhist,se.lesion.hist,'-','color', [.95 .45 .1],'linewidth',2); hold on
    plot(se.nolesion.xhist,se.nolesion.hist,'-','linewidth',2, 'color', [.1 .45 .95]);
    
    % plot the mean lines
    plot([se.nolesion.rmse.mean,se.nolesion.rmse.mean], [0, ymax + yjit],'-','color',[.1 .45 .95] );
    plot([se.lesion.rmse.mean,se.lesion.rmse.mean], [0, ymax + yjit], '-', 'color',[.95 .45 .1]);
    
    % format the plot
    title(sprintf('mean RMSE\nno-lesion: %2.3f | lesion: %2.2f', ...
        se.nolesion.rmse.mean, se.lesion.rmse.mean), 'fontsize', 16);
    ylabel('Probability', 'fontsize', 14);
    xlabel('RMSE', 'fontsize', 14)
    
    % define legend and axes
    legend({'Lesion','No lesion'},'fontsize',14);
    set(gca,'box','off','xtick',[ mindat round(se.xrange(2)/2) round(se.xrange(2)) ], ...
        'ytick', [ 0 ymid ymax ], 'ylim', [ 0 ymax + yjit ], ...
        'xtick', [ mindat middat maxdat ], 'xlim', [ mindat maxdat ], ...
        'tickdir', 'out', 'ticklength', [ 0.025 0 ]);
    
    % increment handle counter
    fi = fi + 1;
    
end

%% bootstraps from strength of evidence

if (strcmp(plt, '2') || strcmp(plt, 'all'))
    
    % extract data from se
    ywo_e = se.s.lesioned_e;
    y_e   = se.s.unlesioned_e;
    woxhis = se.s.lesioned.xbins;
    xhis   = se.s.unlesioned.xbins;
    min_x = se.s.min_x;
    max_x = se.s.max_x;
    
    % drop negative values from y's to fix patch
    y_e(y_e < 0) = 0;
    ywo_e(ywo_e < 0) = 0;
    
    % find maximum value and create a boundary that's higher
    pmax = max(max([y_e(:), ywo_e(:)]));
    bmax = pmax + (pmax * .10);
    
    % plot the shaded distributions
    fig(fi).name = 'virtual_lesion_test_mean_rmse_hist';
    fig(fi).h    = figure('name', fig(fi).name, 'color', 'w');
    fig(fi).type = 'eps';
    set(fig(fi).h, 'Units', 'normalized', 'Position', [0.007 0.55  0.28 0.36]);
    patch([xhis, xhis], y_e(:), [.1 .45 .95], 'FaceColor', [.1 .45 .95], 'EdgeColor', [ 0 0 0 ]); hold on
    patch([woxhis, woxhis], ywo_e(:), [.95 .45 .1], 'FaceColor', [.95 .45 .1], 'EdgeColor', [ 0 0 0 ]);
    
    % format the axes
    set(gca,'tickdir','out', ...
        'box','off', ...
        'ylim',[ 0 bmax ], ...
        'xlim',[ min_x, max_x ], ...
        'ytick',[ 0 pmax/2 pmax ], ...
        'xtick',round(linspace(min_x, max_x, 4)), ...
        'fontsize', 16)
    ylabel('Probability', 'fontsize', 16)
    xlabel('RMSE', 'fontsize', 16)
    title(sprintf('Strength of connection evidence: %2.3f',(se.s.mean)), 'FontSize', 16);
    
    % increment handle counter
    fi = fi + 1;
    
end

%% all VL computed measures

if (strcmp(plt, '3') || strcmp(plt, 'all'))
    
    % create header of plot
    fig(fi).name = 'Size_of_effect_of_the_lesion';
    fig(fi).h    = figure('name',fig(fi).name,'color','w');
    set(fig(fi).h, 'Units', 'normalized', 'Position', [ 0.0035 0.30 0.50 0.50 ]);
    
    % create subplots for each VL measure
    
    % strength of evidence (d-prime)
    subplot(1, 4, 1);
    err = [ se.s.mean, se.s.mean ] + [ -se.s.std, se.s.std ];
    plot(1, se.s.mean, '-o', 'color', [ .95 .45 .1 ], 'linewidth', 2); hold on
    plot([ 1, 1 ], err, '-', 'color', [ .95 .45 .1 ]);
    ylabel('S (s.d.)', 'fontsize', 14);
    set(gca, 'box', 'off', 'xlim', [ 0 2 ], 'xtick', 1, 'XTickLabel', [], ...
        'ylim', [ 0 ceil(se.s.mean + se.s.std) ], ...
        'tickdir', 'out', 'ticklength', [ 0.025 0 ]);
    
    % earth mover's distance
    subplot(1, 4, 2);
    plot(1, se.em.mean, '-o', 'color', [ .95 .45 .1 ], 'linewidth', 2); hold on
    ylabel('Earth mover''s distance (raw scanner units)', 'fontsize', 14);
    set(gca, 'box', 'off', 'xlim', [ 0 2 ], 'xtick', 1, 'XTickLabel', [], ...
        'ylim', [ 0 ceil(se.em.mean) ], ...
        'tickdir', 'out', 'ticklength', [ 0.025 0 ]);
    
    % kl divergence
    subplot(1, 4, 3);
    plot(1, se.kl.mean, '-o', 'color', [ .95 .45 .1 ], 'linewidth', 2); hold on
    ylabel('K-L divergence (bits)', 'fontsize', 14);
    set(gca, 'box', 'off', 'xlim', [ 0 2 ], 'xtick', 1, 'XTickLabel', [], ...
        'ylim', [ 0 max(se.kl.mean)*2 ], ...
        'tickdir', 'out', 'ticklength', [ 0.025 0 ]);
    
    % jeffery's divergence
    subplot(1, 4, 4);
    plot(1, se.j.mean, '-o', 'color', [ .95 .45 .1 ], 'linewidth', 2); hold on
    ylabel('Jeffrey''s divergence (bits)', 'fontsize', 14);
    set(gca, 'box', 'off', 'xlim', [ 0 2 ], 'xtick', 1, 'XTickLabel', [], ...
        'ylim', [ 0 max(se.j.mean)*2 ], ...
        'tickdir', 'out', 'ticklength', [ 0.025 0 ]);
    
end

end

