function [ fh ] = matrix_quick_plot(mat, cmap)
%matrix_quick_plot() returns a fh to a DK matrix plot for fast sanity checks 
%   cmap is the color range - 
%   [0, 4] is what I've used for EMD on test data
%   [-9, -2] is what I've used for Density on test data
%
% Brent McPherson (c), 2017
%

fh = figure('Position', [580 580 1080 680]);
colormap('hot');
imagesc(log(mat));
axis('square'); axis('equal'); axis('tight');
%title('Full Connectome Network');
%xlabel('FS DK Regions');
%ylabel('FS DK Regions');
y = colorbar;
caxis(cmap);
%ylabel(y, 'Log Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

end

