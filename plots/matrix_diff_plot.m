function [ fh, cbar ] = matrix_diff_plot(mat1, mat2)
%matrix_quick_plot() returns a fh to a DK matrix plot for fast sanity checks 
%   cmap is the color range - 
%   [0, 4] is what I've used for EMD on test data
%   [-9, -2] is what I've used for Density on test data
%
% Brent McPherson (c), 2017
%

% calculate the difference, find caxis range
mat = mat1 - mat2;
crng = max(abs(round(minmax(mat(:)'))));

% make the plot
fh = figure('Position', [580 580 1080 680]); 
colormap(redblue);
imagesc(mat);
axis('square'); axis('equal'); axis('tight');
cbar = colorbar;
caxis([-crng crng]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 0]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 0]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 0]);

end

