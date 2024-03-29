function [ fh, mat ] = plotDifferenceMatrix(mat1, mat2, crng)
%plotDifferenceMatrix returns a fh to the difference between mat1 and mat2
% This is used for inspecting differences between the cleaned networks
%
% Brent McPherson (c), 2017
%

% calculate the difference, find caxis range
mat = mat1 - mat2;

if ~exist('crng', 'var') || isempty(crng)
    value = max(abs(round(minmax(mat(:)'))));
    crng = [-value value];
end

% make the plot
fh = figure('Position', [580 580 1080 680]); 
colormap(redblue);
imagesc(mat);
axis('square'); axis('equal'); axis('tight');
colorbar;
caxis(crng);
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);

end

%% accessory fxn of colormap

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

end