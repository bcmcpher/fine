function [ fh ] = plotDifferenceMatrix(mat1, mat2)
%plotDifferenceMatrix() returns a fh to the difference between mat1 and mat2
% This is used for inspecting differences between the cleaned networks
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
colorbar;
caxis([-crng crng]);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);

%% left / right / diagonal lines

full = size(mat, 1) + 0.5;
half = (size(mat, 1) / 2) + 0.5;

line([half half], [0.5 full], 'Color', [0 0 1]);
line([0.5 full], [half half], 'Color', [0 0 1]);
line([full 0.5], [full 0.5], 'Color', [0 0 1]);

end

