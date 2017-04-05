function [ fh ] = plotAdjacencyMatrix(mat, crng)
%plotAdjacencyMatrix() returns a fh to a matrix plot 
%   crng is the color range
%   add option for doing log or not
%   figure out how to add labels w/ fh
%
%   [0, 4] is what I've used for EMD on test data
%   [-2, 7] is what I've used for Density on test data
%
% Brent McPherson (c), 2017
%

fh = figure('Position', [580 580 1080 680]);
colormap('hot');
imagesc(log(mat));
axis('square'); axis('equal'); axis('tight');
colorbar;
caxis(crng);

set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);

%% left / right / diagonal lines

full = size(mat, 1) + 0.5;
half = (size(mat, 1) / 2) + 0.5;

line([half half], [0.5 full], 'Color', [0 0 1]);
line([0.5 full], [half half], 'Color', [0 0 1]);
line([full 0.5], [full 0.5], 'Color', [0 0 1]);

end

