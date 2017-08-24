function [ fig ] = plotVirtualLesion(se)
%plotVirtualLesion recreates the distributions from the LiFE paper for any
% virtual lesion.
%

%% raw RMSE distirbutions
fig(1).name = 'rmse_distributions';
fig(1).h    = figure('name', fig(1).name, 'color', 'w');
fig(1).type = 'eps';

set(fig(1).h,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);

% better relative x-boundary
dat = [ se.nolesion.rmse.all se.lesion.rmse.all ];
jitter = (min(dat) * 0.05);
mindat = round(min(dat) - jitter);
maxdat = round(max(dat) + jitter);
middat = round(maxdat / 2);

% better relative y-boundary
ymax = max([max(se.lesion.hist) max(se.nolesion.hist)]);
ymid = ymax / 2;
yjit = ymax * 0.10;

plot(se.lesion.xhist,se.lesion.hist,'-','color', [.95 .45 .1],'linewidth',2); hold on
plot(se.nolesion.xhist,se.nolesion.hist,'-','linewidth',2, 'color', [.1 .45 .95])
plot([se.nolesion.rmse.mean,se.nolesion.rmse.mean], [0, ymax + yjit],'-','color',[.1 .45 .95] )
plot([se.lesion.rmse.mean,se.lesion.rmse.mean], [0, ymax + yjit], '-', 'color',[.95 .45 .1])

title(sprintf('mean RMSE\nno-lesion %2.3f | lesion %2.2f', se.nolesion.rmse.mean, se.lesion.rmse.mean), 'fontsize', 16)
ylabel('Probability', 'fontsize', 14);
xlabel('RMSE', 'fontsize', 14)

legend({'Lesion','No lesion'},'fontsize',14);

set(gca,'box','off','xtick',[mindat round(se.xrange(2)/2) round(se.xrange(2))],'ytick',[0 ymid ymax],...
    'xlim', [mindat maxdat], 'ylim', [0 ymax + yjit], ...
    'tickdir', 'out', 'ticklength', [0.025 0]);

%% bootstraps from strength of evidence

ywo_e = se.s.lesioned_e;
y_e   = se.s.unlesioned_e;
woxhis = se.s.lesioned.xbins;
xhis   = se.s.unlesioned.xbins;
min_x = se.s.min_x;
max_x = se.s.max_x;

% find maximum of plot
pmax = max(max([y_e(:), ywo_e(:)]));
ran = linspace(0, pmax, 3);

fig(2).name = 'virtual_lesion_test_mean_rmse_hist';
fig(2).h    = figure('name', fig(2).name, 'color', 'w');
fig(2).type = 'eps';
set(fig(2).h, 'Units', 'normalized', 'Position', [0.007 0.55  0.28 0.36]);
%patch([xhis, xhis], y_e(:), [.1 .45 .95], 'FaceColor', [.1 .45 .95], 'EdgeColor', [.1 .45 .95]); hold on
line(xhis, y_e(:, 1)); hold on
line(xhis, y_e(:, 2));
fill(xhis, y_e(:, 2), [.1 .45 .95]);
fill(xhis, y_e(:, 1), 'w');
%patch([woxhis, woxhis], ywo_e(:), [.95 .45 .1], 'FaceColor', [.95 .45 .1], 'EdgeColor', [.95 .45 .1]);
line(woxhis, ywo_e(:, 1));
line(woxhis, ywo_e(:, 2));
fill(woxhis, ywo_e(:, 2), [.95 .45 .1]);
fill(woxhis, ywo_e(:, 1), 'w');

%'ylim',[ran(1) ran(3)], ...
%'ytick', ran, ...
set(gca,'tickdir','out', ...
    'box','off', ...
    'ylim',[0 0.25], ...
    'xlim',[min_x, max_x], ...
    'ytick',[0 0.1 0.2], ...
    'xtick',round(linspace(min_x, max_x, 4)), ...
    'fontsize', 16)
ylabel('Probability', 'fontsize', 16)
xlabel('rmse', 'fontsize', 16)
title(sprintf('Strength of connection evidence %2.3f',(se.s.mean)), 'FontSize', 16);

%% all VL computed measures

fig(3).name = 'Size_of_effect_of_the_lesion';
fig(3).h    = figure('name',fig(3).name,'color','w');
set(fig(3).h,'Units','normalized','Position',[0.007 0.55  0.28 0.36]);
subplot(1,4,1)
plot(1,se.s.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
plot([1,1], [se.s.mean,se.s.mean] + [-se.s.std,se.s.std], '-','color',[.95 .45 .1] )
ylabel('S (s.d.)', 'fontsize',14);
set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.s.mean + se.s.std)], ...
    'tickdir', 'out', 'ticklength', [0.025 0])
subplot(1,4,2)
plot(1,se.em.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
ylabel('Earth mover''s distance (raw scanner units)', 'fontsize',14);
set(gca,'box','off','xlim',[0 2], 'ylim',[0 ceil(se.em.mean)], ...
    'tickdir', 'out', 'ticklength', [0.025 0])
subplot(1,4,3)
plot(1,se.kl.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
ylabel('K-L divergence (bits)', 'fontsize',14);
set(gca,'box','off','xlim',[0 2], 'ylim',[0 max(se.kl.mean)], ...
    'tickdir', 'out', 'ticklength', [0.025 0])
subplot(1,4,4)
plot(1,se.j.mean,'-o','color', [.95 .45 .1],'linewidth',2); hold on
ylabel('Jeffrey''s divergence (bits)', 'fontsize',14);
set(gca,'box','off','xlim',[0 2], 'ylim',[0 max(se.j.mean)], ...
    'tickdir', 'out', 'ticklength', [0.025 0])

end

