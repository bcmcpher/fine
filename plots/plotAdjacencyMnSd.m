function [ fh, ax1, ax2 ] = plotAdjacencyMnSd(mat1, mat2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% check inputs

% if ~exist('names', 'var') || isempty(names)
%     names = [];
% end

% if the matrices aren't the same size, error
if size(mat1) ~= size(mat2)
    error('The inputs have to be the same size');
end

% grab the dimensions
dims = size(mat1, 1);

%% create plots

fh = figure; 

% create the lower left set of axes
ax1 = axes;
imagesc(tril(mat1, 1), 'AlphaData', tril(ones(dims)));
axis('square'); axis('equal'); axis('tight');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', [], 'box', 'off');

% create the upper right set of axes
ax2 = axes;
imagesc(triu(mat2, 1), 'AlphaData', triu(ones(dims)));
axis('square'); axis('equal'); axis('tight');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', [], 'box', 'off');

%% combine axes into 1 figure

% combine into a single figure
linkaxes([ ax1 ax2 ]);

% remove lower axes values
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

% % plot names if passed
% if ~isempty(names)
%     set(gca, 'XTickLabel', names(1:2:dims), 'XTickLabelRotation', 90, 'XTick', 1:2:dims, 'TickDir', 'out');
%     set(gca, 'YTickLabel', names(2:2:dims), 'YTickLabelRotation', 90, 'YTick', 2:2:dims, 'TickDir', 'out', 'box', 'off');
% else
%     set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', [], 'box', 'off');
% end

% set the color maps
colormap(ax1, 'hot');
colormap(ax2, 'parula');

%set([ ax1, ax2 ], 'Position', [ 850 325 800 650 ]);
cb1 = colorbar(ax1, 'Position', [ .05 .11 .0675 .815 ]);
ylabel(cb1, 'Mean');

cb2 = colorbar(ax2, 'Position', [ .88 .11 .0675 .815 ]);
ylabel(cb2, 'SD');

end
