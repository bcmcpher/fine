function [ fh ] = fnPlotBrainNetwork(centers, cl, sz, mat, lc, nnrm, scale)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%
% document / try and break a bit more
% add text labels of the N highest node values?
% lc / nnrm should be the same variable and parsed based on shape?
% be able to set color defaults when weird things are empty
%

%% parse inputs to define options

% grab the number of nodes
nrois = size(centers, 1);

% grab unique, upper diagonal edge indices and number of edges
eind = nchoosek(1:nrois, 2);
edge = size(eind, 1);

% if node color isn't passed, set all to default color
if(~exist('cl', 'var') || isempty(cl))
    cl = repmat([ .17 .17 .34 ], [ nrois 1 ]);
end

% if node scale isn't passed, set all to the same size
if(~exist('sz', 'var') || isempty(sz))
    sz = ones(nrois, 1);
end

% if mat of nrois edge weights isn't passed, set to empty
if(~exist('mat', 'var') || isempty(mat))
    mat = nan(nrois, nrois);
end

% if mat of nrois edge weights isn't passed, set to empty
if(~exist('lc', 'var') || isempty(lc))
    lc = repmat([ 0 0 0 ], [ edge 1 ]);
end

% if nnrm of nrois edge profiles isn't passed, set to empty
if(~exist('nnrm', 'var') || isempty(nnrm))
    ewtp = 'ew';
    nnrm = nan(nrois, nrois, 100);
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
                    
                    % grab the number of nodes
                    nnode = size(y, 1);
                    
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

% clear ii n1 n2 val wdth xs ys zs y xyz tp

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

