function [ fh ] = matrix_quick_plot(mat)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fh = figure('Position', [580 580 1080 680]);
colormap('hot');
imagesc(log(mat));
axis('square'); axis('equal'); axis('tight');
%title('Full Connectome Network');
xlabel('FS DK Regions');
ylabel('FS DK Regions');
y = colorbar;
%ylabel(y, 'Log Number of Streamlines');
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);
line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);

end

