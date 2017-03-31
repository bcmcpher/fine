function [ fh, fgOut ] = render_connection(fe, pconn, label, tract, persp, slice, conn)
% render_connection() renders an ENCODE object with network connection structure  
%   Using a fit FE fiber group and pconn indexes, this fxn renders a
%   connection, path neighborhood, or both with or without cleaning using MBA.
%
% Examples:
% [ fh, fg ] = render_connection(fe, pconn, 2, 'lh-sag', 'fibers', 'path', 0.05);
% [ fh, fg ] = render_connection(fe, pconn, 2, 'axial', 'nzfibs', 'neighborhood', 0.05);
% [ fh, fg ] = render_connection(fe, pconn, 2, 'coronal', 'cln-nzf', 'both', 0.05);
%
% Brent McPherson (c), 2017
%

% load the anatomy from fe
anatomy = niftiRead(fe.path.anatomy);

%% other arguments...

% percent of path neighborhood to display
pnprc = 0.025;

% set minimum length of 10mm
threshold_length = 15;

% generate the blue / red colors for the tracts
colors = {[.1 .25 .65], [.75 .25 .1]};

%% run the analysis

% pull name from pconn based on index
roi1 = pconn{tract}.roi1;
roi2 = pconn{tract}.roi2;
fgName = ['Index: ' num2str(tract) '; ' num2str(roi1) '-to-' num2str(roi2)];

display('Converting streamlines to ACPC space...');
fg = feGet(fe, 'fg acpc');

% create fiber group object for connection
ind1 = getfield(pconn{tract}, label, 'indices');

% if an emtpy fg is requested
if isempty(ind1)
    error('The connection requested contains no streamlines at this parameter.');
end

% create path fiber group
fgOut{1} = fgCreate('name', fgName, 'colorRgb', colors{1}, 'fibers', fg.fibers(ind1));

%% render figure

display(['Rednering ' fgName '...']);

% if I need the path neighborhood
switch conn
    case{'path'}
        
        display('Rendering path...');
        
    case{'neighborhood', 'both'}
        
        display('Rending path neighborhood...');
        
        % pull path neighborhood indices
        ind2 = feGet(fe, 'Path Neighborhood', ind1);
        
        % subset path neighborhood from whole brain fg
        fgOut{2} = fgExtract(fg, ind2, 'keep');
        
        % filter out short path neighborhood streamlines
        c = 1;
        fibers = cell(1, length(fgOut{2}.fibers));
        for ii = 1:length(fgOut{2}.fibers)
            if length(fgOut{2}.fibers{ii}) > threshold_length
                fibers{c} = fgOut{2}.fibers{ii};
                c = c + 1;
            end
        end
        
        fgOut{2}.fibers = fibers;
        clear fibers
        
        % pull a random subset of the path neigborhood to minimize the rendering
        fibs_indx = randsample(1:length(fgOut{2}.fibers), round(length(fgOut{2}.fibers) * pnprc));
        fgOut{2}.fibers = fgOut{2}.fibers(fibs_indx);
        
    otherwise
        error('Invalid connection type requested. Can render either: ''path'', ''neighborhood'', or ''both''.');
        
end

clear fg

%% plot the requested connection / view / fg

for ii = 1:length(persp)
    
    tpersp = persp{ii};
    
    % switch statement to define the plots perspective
    switch tpersp
        case{'lh-sag'}
            viewCoords = [-90, 0];
            indx = [slice(ii) 0 0];
        case{'rh-sag'}
            viewCoords = [90, 0];
            indx = [slice(ii) 0 0];
        case{'axial'}
            viewCoords = [0, 90];
            indx = [0 0 slice(ii)];
        case{'coronal'}
            viewCoords = [180, 0];
            indx = [0 slice(ii) 0];
        otherwise
            error('Invalid slice view requested.');
    end
    
    % plot the brain slice
    fh = mbaDisplayBrainSlice(anatomy, indx, gca); hold on;
    
end

% switch statement for which connections to model
switch conn
    case{'path'}
        [fh, light_h] = mbaDisplayConnectome(fgOut{1}.fibers, fh, colors{1}, 'single');
        delete(light_h);
    case{'neighborhood'}
        [fh, light_h] = mbaDisplayConnectome(fgOut{2}.fibers, fh, colors{2}, 'single');
        delete(light_h);
    case{'both'}
        [fh, light_h1] = mbaDisplayConnectome(fgOut{1}.fibers, fh, colors{1}, 'single');
        delete(light_h1);
        [fh, light_h2] = mbaDisplayConnectome(fgOut{2}.fibers, fh, colors{2}, 'single');
        delete(light_h2);      
end

% plot this way every time
view(viewCoords);
camlight right
lighting phong

end
