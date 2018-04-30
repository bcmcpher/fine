function [ fh, lh, fgConn, fgNeig ] = render_fg(fe, indx, anat, persp, slice, conn)
% render_connection() renders an ENCODE object with network connection structure  
%   Using a fit FE fiber group and pconn indexes, this fxn renders a
%   connection, path neighborhood, or both with or without cleaning using MBA.
%
% Examples:
% [ fh, fg ] = render_connection(fe, indx, {'lh-sag'}, [-1], 'path');
% [ fh, fg ] = render_connection(fe, indx, {'axial'}, [15], 'path');
% [ fh, fg ] = render_connection(fe, indx, {'rh-sag', 'coronal'}, [5, 15], 'neighborhood');
% [ fh, fg ] = render_connection(fe, indx, {'coronal'}, [-15], 'both');
%
% Brent McPherson (c), 2017
%

% load the anatomy from fe
%anat = niftiRead(fe.path.anatomy);

%% other arguments...

% percent of path neighborhood to display
pnprc = 0.025;

% set minimum length of 10mm
threshold_length = 0;

% generate the blue / red colors for the tracts
colors = {[.1 .25 .65], [.75 .25 .1]};

%% run the analysis

% orient fibers
disp('Converting streamlines to ACPC space...');
fg = feGet(fe, 'fg acpc');

% create fiber group object for connection
%ind1 = getfield(pconn{tract}, label, 'indices');

% create path fiber group
for ii = 1:length(indx)
    fgConn{ii} = fgCreate('colorRgb', colors{1}, 'fibers', fg.fibers(indx{ii}));
end

%% render figure

disp('Rendering connection...');

% if I need the path neighborhood
switch conn
    case{'path'}
        
        disp('Rendering path...');
        
    case{'neighborhood', 'both'}
        
        disp('Rending path neighborhood...');
        
        % pull path neighborhood indices
        ind2 = feGet(fe, 'Path Neighborhood', ind1);
        
        % subset path neighborhood from whole brain fg
        fgNeig{1} = fgExtract(fg, ind2, 'keep');
        
        % filter out short path neighborhood streamlines
        c = 1;
        fibers = cell(1, length(fgNeigh{1}.fibers));
        for ii = 1:length(fgNeigh{1}.fibers)
            if length(fgNeigh{1}.fibers{ii}) > threshold_length
                fibers{c} = fgNeigh{1}.fibers{ii};
                c = c + 1;
            end
        end
        
        fgNeigh{1}.fibers = fibers;
        clear fibers
        
        % pull a random subset of the path neigborhood to minimize the rendering
        fibs_indx = randsample(1:length(fgConn{2}.fibers), round(length(fgConn{2}.fibers) * pnprc));
        fgConn{2}.fibers = fgConn{2}.fibers(fibs_indx);
        
    otherwise
        error('Invalid connection type requested. Can render either: ''path'', ''neighborhood'', or ''both''.');
        
end

clear fg

%% plot the requested connection / view / fg

fh = figure;

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
    fh = mbaDisplayBrainSlice(anat, indx, gca); hold on;
    
end

% switch statement for which connections to model
switch conn
    case{'path'}
        for ii = 1:length(indx)
            [fh, light_h] = mbaDisplayConnectome(fgConn{ii}.fibers, fh, colors{1}, 'single', 'fiberRadius', 0.15);
            delete(light_h);
        end
    case{'neighborhood'}
        [fh, light_h] = mbaDisplayConnectome(fgConn{2}.fibers, fh, colors{2}, 'single', 'fiberRadius', 0.15);
        delete(light_h);
    case{'both'}
        [fh, light_h1] = mbaDisplayConnectome(fgConn{1}.fibers, fh, colors{1}, 'single', 'fiberRadius', 0.15);
        delete(light_h1);
        [fh, light_h2] = mbaDisplayConnectome(fgConn{2}.fibers, fh, colors{2}, 'single', 'fiberRadius', 0.15);
        delete(light_h2);      
end

% plot this way every time
view(viewCoords);
lh = camlight('right');
lighting phong

end


% % combine all the nodes of all the streamlines
% fibs = [];
% for jj = 1:length(fibers)
%     fibs = [ fibs fibers{jj} ];
% end
% 
% %  find min / max [x,y,z]
% axCoords = minmax(fibs);
% 
% %  compute 5% padding of [x,y,z]
% axPadding = abs(axCoords(:,1) - axCoords(:,2)) * 0.10;
% 
% % add padding to [x,y,z] dimensions
% axOut = [ round(axCoords(:, 1) + 4*axPadding), round(axCoords(:, 2) - 4*axPadding) ];
% 
% % decide the slice
% indx=[ 0 0 axOut(3,1) ];
% 
% % apply default limits to axis
% axOut2 = [round(axCoords(:,1) - axPadding), round(axCoords(:,2) + axPadding)];
% axis([reshape(axOut2', 1, 6)]);
