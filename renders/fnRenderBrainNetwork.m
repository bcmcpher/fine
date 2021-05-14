function [ fh ] = fnRenderBrainNetwork(netw, nodec, nodez, edgew, thr, edgep, dtype, scale, noaxis)
%[ fh ] = fnRenderBrainNetwork(netw, nodec, nodez, edgew, thr, edgep, dtype, scale, noaxis);
%     This function renders the netw object as a ball-stick brain shaped 
% network summary of the data. You can click on the node in the plots to 
% display the label if they're stored during network creation. 
% 
% INPUTS:
%     netw    - network object
%     nodec   - N x 3 matrix of RGB colors for each node (optional)
%     nodez   - N x 1 vector to denote the scaling of each node
%               (lager value == bigger node) (optional)
%     edgew   - a string requesting the edge weight to render between nodes
%     thr     - the percentile threshold to apply to edges. Between 0 and
%               1 (0 removes no edges, 1 removes all edges)
%     edgep   - edge property; values to display along the edge (optional)
%               M is the number of possible edges. It may be either:
%                 - a M x 3 vector of RGB values of edge assignments
%                 - a string of 'tprof', 'profile', or 'profiles'
%     dtype   - a string requesting the tract profile to render as each edge
%     scale   - a heuristic for the relative scale of the node size and
%               line thickness. default of 5. Change by single digits.
%               (optional)
%     noaxis  - clear background color and axes for a cleaner view to save
%
% OUTPUTS:
%     fh - a figure handle for the plot
%               
% TODO: 
% - clear up variable names / logic
% - better name rendering options?
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

%% parse inputs for basic definition of nodes / edges

% pull centers and labels from network object
nodep = cell2mat(cellfun(@(x) x.center.acpc, netw.nodes, 'UniformOutput', false));
labels = cellfun(@(x) x.name, netw.nodes);

% grab the number of nodes
nnodes = size(nodep, 1);

% grab unique, upper diagonal edge indices and number of edges
pairs = netw.parc.pairs;
nedges = size(netw.edges, 1);

% if node color isn't passed, set all to default color
if(~exist('nodec', 'var') || isempty(nodec))
    nodec = repmat([ .67 .67 .67 ], [ nnodes 1 ]);
end

% if node scale isn't passed, set all to the same size
if(~exist('nodez', 'var') || isempty(nodez))
    nodez = ones(nnodes, 1);
end

% if scale is not passed, set to default
if(~exist('scale', 'var') || isempty(scale))
    scale = 5;
end

% if scale is not passed, set to default
if(~exist('noaxis', 'var') || isempty(noaxis))
    noaxis = true;
end

% if the matrix field isn't made, don't bother trying
if ~isfield(netw.edges{1}, 'matrix')
    error('Network edge weights have not been summarized yet. Run: netw = fnComputeMatrixEdges(netw);');
else
    % if the mat field is not passed, keep empty, don't draw edges
    if(~exist('edgew', 'var') || isempty(edgew))
        edgew = []; % if it's empty, don't draw edges
        draw_edges = false;
    else
        if ~isfield(netw.edges{1}.matrix, edgew)
            warning('The requested field is not a valid edge weight. Setting to ''density''.');
            edgew = 'density';
        end
        draw_edges = true;
    end
end

% set threshold to the 50 percentile if they don't request something else
if(~exist('thr', 'var') || isempty(thr))
    thr = 0.50;
end

% scale threshold to be between 0 and 1
if (thr > 1) || (thr < 0)
    warning('Threshold must be set as a positive value <= 1. Resetting to 0.50');
    thr = 0.50;
end

%% scale the node / edge inputs

% compute normalized node size and scale for display
nodez = (nodez / max(nodez)) * scale;

% if edges are to be drawn
if draw_edges
    
    % pull the requested edge weight
    edgev = cellfun(@(x) x.matrix.(edgew), netw.edges);
    
    % compute a threshold on the matrix edge
    sthr = thr * 100; % scale for percentile to work right
    edgev(edgev < prctile(edgev, sthr)) = 0;
    
    % subset pairs and ewgh to only non-zero values so edge objects are
    % sparse and only contain values getting plotted
    ikeep = edgev ~= 0;
    pairs = pairs(ikeep, :);
    edgev = edgev(ikeep);
    
    % pull the edge indices and get the count of kept edges
    edgek = find(ikeep); 
    nkeep = sum(ikeep);
    
    % normalize edge weights so the line thickness is scaled usefully
    edgev = edgev / max(edgev);
    
else
    
    % just set edge vector as empty
    edgev = [];
    nkeep = [];
    
end

%% deal with edges

% if the edge property (color/profiles) are empty, fill in edge color as black
if(~exist('edgep', 'var') || isempty(edgep))
    
    edgec = repmat([ 0 0 0 ], [ nkeep, 1 ]); % set all edges black
    draw_profiles = false;
    
else % otherwise
    
    % if dtype for profile isn't passed, don't compute profiles
    if(~exist('dtype', 'var') || isempty(dtype))
        if ischar(edgep) && any(strcmp(edgep, {'tprof','profile','profiles'}))
            error('A ''dtype'' label must also be passed if profiles are requested.');
        end
        dtype = [];
        draw_profiles = false;
    else
        draw_profiles = true;
        cmap = parula(100); % define the colormap for the profiles
    end
    
    % if edgep is an nedges x 3 matrix, assume RGB colors
    if size(edgep, 1) == nedges && size(edgep, 2) == 3
        
        edgec = edgep; % use the passed matrix to color the edges
        edgec = edgec(ikeep, :); % subset the passed colors for each edge down to what is plotted
        draw_profiles = false; % override the old setting
        if ~isempty(dtype)
            warning('Based on the arguments passed, the requested dtype ''%s'' is not used in this plot.', dtype);
        end
        
    % otherwise, if edgep is a character string and draw_profiles is expected
    elseif (ischar(edgep) && any(strcmp(edgep, {'tprof','profile','profiles'}))) && draw_profiles 
        
        % if the requested dtype profiles aren't stored, error
        if ~isfield(netw.edges{1}.profile, dtype)
            error('The requested profile ''%s'' is not stored.', dtype);
            
        else
            
            % otherwise build the profile matrix based on input size
            nnode = size(netw.edges{1}.profile.(dtype), 1);
            trpof = nan(nkeep, nnode); % preallocate profiles
            for edge = 1:nkeep
                trpof(edge,:) = netw.edges{edgek(edge)}.profile.(dtype)';
            end
            %ttube = cell(nkeep, 1); preallocate cell array to catch profile tubes
        end
    else % otherwise the combination of variables passed is too confusing
        error('''edgep'' does not have a valid operation for value ''%s''.', edgep);
    end
end

clear edge ikeep

%% define the resolution of the rendered surfaces

% create cannonical sphere - sets the resolution of a node's point
NumSphFaces = 15;
[ SX, SY, SZ ] = sphere(NumSphFaces);

% tube resolution
subdivs = 10;

%% define the plot limits

% grab the image dimension for min / max
bbox = [ 1 1 1; netw.parc.dsize ]';

% convert the image dimensions to acpc space - transposed for easier parsing by plot axis
plim = mrAnatXformCoords(netw.parc.xform.img2acpc, bbox)';

% set the limits of the image frame based on the parc image
% sort to ensure neg/pos order for plotting (RAS vs LAS flips x, etc.)
xlim = sort(plim(1, :)); %[ -80 80 ];
ylim = sort(plim(2, :)); %[ -100 80 ];
zlim = sort(plim(3, :)); %[ -60 80 ];

clear bbox plim

%% create the nodes w/ requested parameters

% create and catch the node surfaces w/ the option to toggle the labels
nsurf = cell(nnodes, 1);

% initialize the plot
fh = figure; hold on;

% plot the spheres - for every node in nnodes
for node = 1:nnodes
    
    % scale / align the canonical sphere to the node center, scale size by input vector
    nsurf{node} = surf(SX*nodez(node) + nodep(node, 1), SY*nodez(node) + nodep(node, 2), SZ*nodez(node) + nodep(node, 3));
    
    % change color and remove edges
    nsurf{node}.FaceColor = nodec(node, :);
    nsurf{node}.EdgeAlpha = 0;
        
    % print passed label, but set it to be invisible by default
    tlab = text(nodep(node, 1), nodep(node, 2), nodep(node, 3), [ '    ' labels{node} ], ...
                'HorizontalAlignment', 'left', 'FontSize', scale*2, 'Visible', 'off');
    
    % create the click option to toggle the label visible / invisible
    nsurf{node}.ButtonDownFcn = {@toggleText, tlab};
    
    % set reflection parameters
    material metal
    
end

clear node tlab SX SY SZ NumSphFaces

%% define how edges will be drawn

% if mat is not entirely empty, plot edges
if ~all(isnan(edgev(:))) && draw_edges
    
    % plot the lines
    for edge = 1:nkeep
        
        % grab the indices / edge value
        n1 = pairs(edge, 1);
        n2 = pairs(edge, 2);
        val = edgev(edge);
        
        % if there's an edge above threshold
        if val > 0
            
            % create node centers / line endpoints
            xs = [ nodep(n1, 1) nodep(n2, 1) ];
            ys = [ nodep(n1, 2) nodep(n2, 2) ];
            zs = [ nodep(n1, 3) nodep(n2, 3) ];
            
            if ~draw_profiles
                
                % compute relative line thickness from edge weight
                wdth = val * scale;
                
                % plot the line w/ edge thickness / color
                plot3(xs, ys, zs, 'color', edgec(edge, :), 'LineWidth', wdth);
                
            else
                
                % pull just the relevant profile
                y = trpof(edge, :);
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
                r = val * scale;
                
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
                
                % If color has a row for each fiber then each fiber node will have its own color
                fcolor = reshape(color, [ N 1 3 ]);
                C = repmat(fcolor,[ 1 subdivs 1 ]);
                
                % render tube
                surf(X, Y, Z, C, 'EdgeColor', 'none');
                % it started drawing yellow lines along the tube...
                
            end
            
        end
        
    end
    
end

clear edge cmap color fcolor node n1 n2 theta val wdth C N X Y Z b jj n r s xs ys zs y xyz tp

% set axes
set(gca, 'xlim', xlim, 'ylim', ylim, 'zlim', zlim);
axis equal

% choose view
view(0, 90);

% fix lighting
camlight('right');
lighting phong

% hide axes if requested
if noaxis
    set(gca, 'color', 'none', 'xtick', [], 'ytick', [], 'ztick', [], 'visible', 'off');
end

hold off

end

% toggle the text on and off
function toggleText(nsurf, event, lab)

% nsurf is nsurf{node} that is clicked, has to be passed here
% event is a Hit that is the click event - no need to interact with it
% lab is the label that is stored within nsurf that is set to invisible by default

% print an obnoxious statement if the node is somehow invisible so it stops
% complaining that I don't use nsurf in this fxn
if strcmp(nsurf.Visible, 'off')
    disp('How did you even see this to select it?');
end

% pull visibility from stored label
vis = get(lab, 'Visible');

% if it's on, toggle off - else, toggle on
if strcmp(vis, 'on')
   disp([ 'Hiding label for node nearest: ' num2str(event.IntersectionPoint) ]);
   set(lab, 'Visible', 'off');
else
   disp([ 'Displaying label for node nearest: ' num2str(event.IntersectionPoint,2) ]);
   set(lab, 'Visible', 'on');
end

end
