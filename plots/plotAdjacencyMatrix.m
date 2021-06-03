function [ fh ] = plotAdjacencyMatrix(mat, names)
%plotAdjacencyMatrix() returns a fh to a matrix plot 
%   crng is the color range
%   add option for doing log or not
%   figure out how to add labels w/ fh
%   indices to sort edges
%   lines that divide modules 
%
%   [0, 4] is what I've used for EMD on test data
%   [-2, 7] is what I've used for Density on test data
%
% Brent McPherson (c), 2017
%

if ~exist('names', 'var') || isempty(names)
    names = [];
end

% if ~exist('lines', 'var') || isempty(lines)
%   lines = 'noshow';
% end

% if ~exist('crng', 'var') || isempty(crng)
%     if max(mat(:)) > 0
%         crng = [ 0 max(mat(:)) ];
%     else
%         crng = [ min(mat(:)) 0 ];
%     end
% end

fh = figure('Position', [580 580 1080 680]);
colormap('hot');
imagesc(mat);
axis('square'); axis('equal'); axis('tight');
colorbar;
%caxis(crng);

if ~isempty(names)
    set(gca, 'XTickLabel', names, 'XTickLabelRotation', 45, 'XTick', 1:size(mat,1), 'TickDir', 'out', 'box', 'off');
    set(gca, 'YTickLabel', '', 'YTick', []);
else
    set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', [], 'box', 'off');
end

% %% left / right / diagonal lines
% 
% if ischar(lines)
%     
%     switch lines
%         
%         case {'equal'}
%             
%             full = size(mat, 1) + 0.5;
%             half = (size(mat, 1) / 2) + 0.5;
%             
%             % these are off a bit
%             line([half half], [0.5 full], 'Color', [1 1 1], 'Linewidth', 0.25);
%             line([0.5 full], [half half], 'Color', [1 1 1], 'Linewidth', 0.25);
%             line([full 0.5], [full 0.5], 'Color', [1 1 1], 'Linewidth', 0.25);
%             
%         case {'noshow'}
%             
%         otherwise
%             
%             warning('Unknown line option requested. No lines drawn.');
%             
%     end
%     
% else
% 
%     % assume this is plotting module boundaries
%     nummod = max(lines);
%     N = size(mat, 1);
%     mod_boundaries = find(abs(diff(lines))>0);
%     for i=1:nummod-1
%         ll = line([mod_boundaries(i)+0.5 mod_boundaries(i)+0.5],[0.5 N+0.5]); set(ll,'Color',[1 1 1],'LineWidth',0.25);
%         ll = line([0.5 N+0.5],[mod_boundaries(i)+0.5 mod_boundaries(i)+0.5]); set(ll,'Color',[1 1 1],'LineWidth',0.25);
%     end
% 
% end

end
