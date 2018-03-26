function [ fh ] = fnPlotModulePoints(mats, clr)
%fnPlotModulePoints() plots each unique module of each subjects as a point
%
%   Hopefully this shows group differences...
%

% pull the number of modules from input
nmod = size(mats, 1);

% unique group should match the # of colors passed
nsubj = size(mats, 3);

% pull maximum range of axes
amax = ceil(max(mats(:)));
amid = round(amax / 2);

% create diagonal indices from input size
diag = [ 1:nmod; 1:nmod ]';

% create the upper diagonal combinations
updg = nchoosek(1:nmod, 2);

% create the index order to plot all the points
indx = [ diag; updg ];

% total number of data points
npts = size(indx, 1);

% create module labels
labs = cell(npts, 1);
for lab = 1:npts
    labs{lab} = [ num2str(indx(lab, 1)) '-' num2str(indx(lab, 2)) ];
end

% define subject jitter
jitter = linspace(-0.20, 0.20, length(1:nsubj));

% plot the data
fh = figure('Position', [ 150 575 1675 625 ]); hold on;

% for every subject
for subj = 1:nsubj
    
    % pull subject data and color from input
    mat = mats(:, :, subj);
    color = clr{subj};
    
    % for every module
    for mod = 1:size(indx, 1)
        
        % the data value to plot
        pt = mat(indx(mod, 1), indx(mod, 2));
        
        % plot the data point
        plot(mod + jitter(subj), pt, 'o', 'MarkerFaceColor', color, ...
            'MarkerEdgeColor', 'black', 'LineWidth', 0.25, 'MarkerSize', 8);
        
    end
    
end

% format the plot axes
set(gca, 'YTick', [ 0 amid amax ], 'YTickLabel', [ 0 amid amax ]);
set(gca, 'XTick', 1:npts, 'XTickLabel', labs, 'XLim', [ 0 npts + 1 ], 'TickLength', [ 0 0 ]);

end

