function [ fh ] = fnPlotBrainNetwork(centers, cl, sz, mat, eprp, label, scale)
%fh = fnPlotBrainNetwork(centers, cl, sz, mat, eprp, label, scale); 
%takes output from the FiNE tools and creates a ball and stick brain-shaped
%network summary of the data.
%
% It is important that the order of the nodes is the same for every input
% that is passed. The simplest way to achieve this is to NOT SORT THE OUTPUTS. 
% 
% INPUTS:
%     centers - N x 3 matrix of AC-PC coordinates representing the center
%               of each cortical node (N is the number of nodes to plot)
%     cl      - N x 3 matrix of RGB colors for each node (optional)
%     sz      - N x 1 vector to denote the scaling of each node
%               (lager value == bigger node) (optional)
%     mat     - N x N connectivity matrix. If a weighted matrix is passed,
%               the values weight the line thickness (optional)
%     eprp    - edge property; values to display along the edge (optional)
%               M is the number of possible edges. It may be either:
%                 - a M x 3 vector of RGB values of edge assignments
%                 - a N x N x nnode matrix of tract profiles
%     label   - N x 1 cell array of text labels for each node. The highest
%               weighted 5% are labeled based on sz. If sz is empty, the 
%               first 5% are labeled (optional)
%     scale   - a heuristic for the relative scale of the node size and
%               line thickness. default of 5. Change by single digits.
%               (optional)
%
% OUTPUTS:
%     fh - a figure handle for the plot
%               
% TODO: 
% - be able to set color defaults when weird things are empty(?)
%   might already be handled OK
%   
% EXAMPLES:
%
% % just plot the nodes
% fh = fnPlotBrainNetwork(centers);
%
% % nodes with color
% fh = fnPlotBrainNetwork(centers, cl);
%
% % nodes with color and scaled
% fh = fnPlotBrainNetwork(centers, cl, sz);
%
% % nodes with color, scale, and edges
% fh = fnPlotBrainNetwork(centers, cl, sz, mat);
%
% % same plot as above w/ colored edges
% fh = fnPlotBrainNetwork(centers, cl, sz, mat, lc);
%
% % same plot as above w/ profiles along edges
% fh = fnPlotBrainNetwork(centers, cl, sz, mat, nnrm);
% 
% % same plots w/ a slightly larger scale of the nodes
% fh = fnPlotBrainNetwork(centers, cl, sz, mat, lc, 7);
% fh = fnPlotBrainNetwork(centers, cl, sz, mat, nnrm, 7);
%
% % add text labels to the nodes
% fh = fnPlotBrainNetwork(centers, cl, sz, [], [], label);
%
% Brent McPherson (c), 2018 - Indiana University 
%

%% parse inputs to define options

% grab the number of nodes
nrois = size(centers, 1);

% grab unique, upper diagonal edge indices and number of edges
eind = nchoosek(1:nrois, 2);
edge = size(eind, 1);

% if node color isn't passed, set all to default color
if(~exist('cl', 'var') || isempty(cl))
    cl = repmat([ .67 .67 .67 ], [ nrois 1 ]);
end

% if node scale isn't passed, set all to the same size
if(~exist('sz', 'var') || isempty(sz))
    sz = ones(nrois, 1);
end

% if mat of nrois edge weights isn't passed, set to empty
if(~exist('mat', 'var') || isempty(mat))
    mat = nan(nrois, nrois);
end

% if the edge color / profiles are empty, it's empty, otherwise
if(~exist('eprp', 'var') || isempty(eprp))
    
    % set to empty
    eprp = [];
    
else
    
    % if it's a edge x 3 matrix of RGB colors
    if size(eprp, 1) == edge && size(eprp, 2) == 3
        lc = eprp;
    end
    
    % if it's a nrois x nrois x nnode matrix w/ tract profile data
    if size(eprp, 1) == nrois && size(eprp, 2) == nrois && size(eprp, 3) > 1
        nnrm = eprp;
        nnode = size(eprp, 3);
    end
    
end

% if edge community colors aren't passed, edges are black
if(~exist('lc', 'var') || isempty(lc))
    lc = repmat([ 0 0 0 ], [ edge 1 ]);
end

% if nnrm of nrois edge profiles isn't passed, set to empty
if(~exist('nnrm', 'var') || isempty(nnrm))
    ewtp = 'ew';
    nnrm = [];
end

% if nnrm of nrois edge profiles isn't passed, set to empty
if(~exist('label', 'var') || isempty(label))
    label = {};
end

% if nnrm exists and mat doesn't exist, throw a warning saying profs aren't plotted
if(~isempty(nnrm) && all(isnan(mat(:))))
    warning('The ''mat'' argument must be passed to plot profile weighted edges.');
end

% if scale is not passed, set to default
if(~exist('scale', 'var') || isempty(scale))
    scale = 5;
end

%% fixed plotting parameters

% create cannonical sphere - sets the resolution of a node's point
NumSphFaces = 15;
[ SX, SY, SZ ] = sphere(NumSphFaces);

% hard set the image bounds - take from volume?
xlim = [ -80 80 ];
ylim = [ -100 80 ];
zlim = [ -60 80 ];

% tube resolution
subdivs = 10;

% color map for potential tract profiles
cmap = hot(100);

%% scale the node / edge inputs

% compute normalized node size and scale for display
% sz = sz * 500; % for scatter3
sz = (sz / max(sz)) * scale;

% sort node size to find top ranked nodes for potential labels
[ ~, si ] = sort(sz, 'descend');

% subset to the top 5% of node degree
nlab = si(1:round(.05 * nrois));

% compute a normalized line thickness

% pull the edge weights
ewgh = nan(edge, 1);
for ii = 1:edge
    ewgh(ii) = mat(eind(ii, 1), eind(ii, 2));
end

% drop if less than zero
%ewgh = ewgh(ewgh > 0);

% grab number to normalize with
ewmx = max(ewgh(:));

%% create the nodes w/ requested parameters

% plot
fh = figure; hold on;

% plot the spheres

% for every K node
for K = 1:nrois
    
    % scale / align the canonical sphere to the node center, scale size by input vector
    tsf = surf(SX*sz(K) + centers(K, 1), SY*sz(K) + centers(K, 2), SZ*sz(K) + centers(K, 3));
    
    % change color and remove edges
    tsf.FaceColor = cl(K, :);
    tsf.EdgeAlpha = 0;
    
    % set reflection parameters
    material metal
    
end

clear K tsf

%% define how edges will be drawn

% edge weight or tract profile
%ewtp = 'tp';

% if ewtp is still empty, assume profiles were passed
if(~exist('ewtp', 'var') || isempty(ewtp))
    ewtp = 'tp';
end

% if mat is not entirely empty, plot edges
if ~all(isnan(mat(:)))
    
    % plot the lines
    for ii = 1:edge
        
        % grab the indices / edge value
        n1 = eind(ii, 1);
        n2 = eind(ii, 2);
        val = mat(n1, n2);
        
        % if there's an edge above threshold
        if val > 0
            
            % create node centers / line endpoints
            xs = [ centers(n1, 1) centers(n2, 1) ];
            ys = [ centers(n1, 2) centers(n2, 2) ];
            zs = [ centers(n1, 3) centers(n2, 3) ];
            
            switch ewtp
                
                case 'ew'
                    
                    % compute relative line thickness from edge weight
                    wdth = (val / ewmx) * 7;
                    
                    % plot the line w/ edge thickness / color
                    plot3(xs, ys, zs, 'color', lc(ii, :), 'LineWidth', wdth);
                    
                case 'tp'
                    
                    % pull just the relevant profile
                    y = squeeze(nnrm(n1, n2, :));
                    y(isinf(y)) = 0;
                    y(isnan(y)) = 0;
                    
                    % skip empty profiles
                    if all(y == 0)
                        continue
                    end
                    
                    % grab each tp node color mapped value
                    nrmy = floor((y / max(y))*100);
                    
                    % skip any NaN values
                    if any(isnan(nrmy)) || any(nrmy > 100)
                        continue
                    end
                    
                    % for each edge, create a uniquely ranged color map?
                    color = cmap(nrmy, :);
                    
                    % create the spacing of points between
                    xyz = [ linspace(xs(1), xs(2), nnode).' linspace(ys(1), ys(2), nnode).' linspace(zs(1), zs(2), nnode).' ];
                    
                    % radius scale with edge weight
                    r = (val / ewmx) * scale;
                    
                    % preallocate mesh coords for tube
                    N = size(xyz, 1);
                    X = zeros(N, subdivs);
                    Y = zeros(N, subdivs);
                    Z = zeros(N, subdivs);
                    theta = 0:(2 * pi / (subdivs - 1)):(2 * pi);
                    
                    % external fxn in MBA and AFQ
                    % creates normal / binormal points for rendering cylinder
                    [ ~, n, b] = frame(xyz(:, 1), xyz(:, 2), xyz(:, 3), randn(1, 3));
                    
                    % set the radius of the tube
                    r = r * ones(N, 1);
                    
                    % Build a mesh for the fiber
                    for jj = 1:N
                        X(jj, :) = xyz(jj,1) + r(jj)*(n(jj,1)*cos(theta) + b(jj,1)*sin(theta));
                        Y(jj, :) = xyz(jj,2) + r(jj)*(n(jj,2)*cos(theta) + b(jj,2)*sin(theta));
                        Z(jj, :) = xyz(jj,3) + r(jj)*(n(jj,3)*cos(theta) + b(jj,3)*sin(theta));
                    end
                    
                    % If color has a row for each fiber than each fiber node will have its own color
                    fcolor = reshape(color, [ N 1 3 ]);
                    C = repmat(fcolor,[ 1 subdivs 1 ]);
                    
                    % render tube
                    surf(X, Y, Z, C, 'EdgeColor', 'interp');
                    
                otherwise
                    error('Not a valid edge request.');
                    
            end
            
        end
        
    end

end

% add the labels if they exist
if ~isempty(label)
    
    for K = 1:nrois
        
        % if it's in the top 5%
        if any(nlab == K)
            
            % print passed label
            text(centers(K, 1), centers(K, 2), centers(K, 3), [ '    ' label{K} ], ...
                 'HorizontalAlignment', 'left', 'FontSize', 10);
            
        end
        
    end
end

clear ii K n1 n2 val wdth xs ys zs y xyz tp

% set axes
set(gca, 'xlim', xlim, 'ylim', ylim, 'zlim', zlim);
axis equal

% choose view
view(0, 90);

% fix lighting
camlight('right');
lighting phong

hold off

end

