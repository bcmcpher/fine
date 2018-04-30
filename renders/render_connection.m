function [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, label, tract, anat, color, persp, slice, conn, crop)
%[ fh, lh, fgOut, pnOut, xyz ] = render_connection() renders any number of  
% specified connections w/ or w/o brain slices for anatomical
% visualizations.
%
% this function requires: parpool
%
% INPUTS:
%     fe      - an evaluated fe structure
%     pconn   - paired connection object storing streamline edge assignments
%     label   - string indicating the fiber groups for which to create tract profiles
%               either:
%                   'all' for all assigned streamlines or
%                   'nzw' for non-zero weighted fibers returned by LiFE
%                   'zwr' for all zero weighted removed streamlines
%     tract   - the index of connections within pconn to render
%     anat    - a loaded nifti object for displaying background brain slices
%     color   - a cell array of length(tract) that contains RGB
%               designations for each color. The default assigns blue to
%               each connection.
%     persp   - a cell array of the of slice orientations to display. 
%               The final orientation is the view that is returned.
%               options are:
%                   'lh-sag'  - left sagittal view
%                   'rh-sag'  - right sagittal view
%                   'axial'   - top-down axial view
%                   'coronal' - coronal view
%     slice   - an array of length(persp) indicating the mm index from the
%               AC-PC point that the slice will be displayed
%     conn    - what streamlines from the requested connections are displayed
%               either: 
%                   'tract' or 'path'
%                   'neighborhood'
%                   'both'
%     crop    - by default, crop the view to the requested streamlines and slices
%             
% OUTPUTS:
%     fh    - figure handle
%     lh    - light handle
%     fgOut - a cell array of the displayed fiber groups
%     pnOut - a cell array of the fiber group path neighborhoods
%             these are randomly pruned to a certain percent to ease rendering
%     xyz   - the initial axis() call to crop view to streamlines / slices
%
% TODO:
% - add tract profile tube
%
% Examples:
% 
% color = {[ .2 .8 .2]};
% views = {'axial', 'lh-sag'};
% indx = [ 5 -35 ];
% 
% [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, 'nzw', 13);
% [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, 'nzw', 13, nii, {[.2, .8, .2]});
% [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, 'nzw', 13, nii, ...
%                                                   {[.2, .8, .2]}, {'axial', 'lh-sag'}, indx);
% [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, 'nzw', 13, nii, color, views, indx, 'neighborhood');
% [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, 'nzw', 13, nii, color, views, indx, 'both');
% [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, 'nzw', 13, nii, color, views, indx, 'tract', 0);
% [ fh, lh, fgOut, pnOut, xyz ] = render_connection(fe, pconn, 'nzw', 13, nii, color, views, indx, 'tract', 1);
%
% Brent McPherson (c), 2018, Indiana University
%

%% other arguments...?

% percent of path neighborhood to display
pnprc = 0.025;

% set minimum length of 10mm
threshold_length = 15;

% generate the blue / red colors for the tracts
colors = {[.1 .25 .65], [.75 .25 .1]};

%% parse the arguments

% assume slices are plot
do_slice = 1;

% try to load the anatomy from fe if anatomy is empty
if(~exist('anat', 'var') || isempty(anat))
    try
        disp('Loaded anatomy file is not passed. Trying to load it from fe path...');
        anat = niftiRead(fe.path.anatomy);
    catch
        warning('''anat'' argument is not found and cannot be loaded. No brain slices will be displayed.');
        do_slice = 0;
    end
end

% if color is empty, set all to blue
% all path neighborhoods are the same red
if(~exist('color', 'var') || isempty(color))
    [ color{1:length(tract)} ] = deal(colors{1});
end
    
% make sure slice is optional
if(~exist('persp', 'var') || isempty(persp))
    persp = 'axial';
    do_slice = 0;
end

% make sure slice is optional
if(~exist('slice', 'var') || isempty(slice))
    slice = [];
    do_slice = 0;
end

% conn defaults to tract
if(~exist('conn', 'var') || isempty(conn))
    conn = 'tract';
end

% crop defaults to on
if(~exist('crop', 'var') || isempty(crop))
    crop = 1;
end

%% run the analysis

disp('Converting streamlines to AC-PC space...');
fg = feGet(fe, 'fg acpc');

% preallocate tracts and path neighborhoods
fgOut = cell(length(tract), 1);
pnOut = cell(length(tract), 1);

% for every tract, compute both
for tct = 1:length(tract)
    
    % pull name from pconn based on index
    roi1 = pconn{tract(tct)}.roi1;
    roi2 = pconn{tract(tct)}.roi2;
    fgName = ['Index: ' num2str(tract(tct)) '; ' num2str(roi1) '-to-' num2str(roi2)];
    
    % create fiber group object for connection
    ind1 = getfield(pconn{tract(tct)}, label, 'indices');
    
    % if an emtpy fg is requested
    if isempty(ind1)
        warning(['Requested connection ' num2str(tract(tct)) ' contains no streamlines.' ]);
    end
    
    disp([ 'Extracting path ' num2str(tct) '...' ]);
    
    % create path fiber group
    fgOut{tct} = fgCreate('name', fgName, 'colorRgb', color{tct}, 'fibers', fg.fibers(ind1));
    
    disp([ 'Extracting path neighborhood ' num2str(tct) '...' ]);
    
    % pull path neighborhood indices
    ind2 = feGet(fe, 'Path Neighborhood', ind1);
    
    % subset path neighborhood from whole brain fg
    pnOut{tct} = fgExtract(fg, ind2, 'keep');
    
    % filter out short path neighborhood streamlines
    c = 1;
    fibers = cell(1, length(pnOut{tct}.fibers));
    for ii = 1:length(pnOut{tct}.fibers)
        if length(pnOut{tct}.fibers{ii}) > threshold_length
            fibers{c} = pnOut{tct}.fibers{ii};
            c = c + 1;
        end
    end
    
    % add the the output fg array
    pnOut{tct}.fibers = fibers;
    clear fibers
    
    % pull a random subset of the path neigborhood to minimize the rendering
    fibs_indx = randsample(1:length(pnOut{tct}.fibers), round(length(pnOut{tct}.fibers) * pnprc));
    pnOut{tct}.fibers = pnOut{tct}.fibers(fibs_indx);
    
end

clear fg

%% plot the requested connection / view / fg

fh = figure;

xbnd = [];

if do_slice

    disp(['Displaying ' num2str(length(tract)) ' brain slices...']);
    
    for ii = 1:length(persp)
        
        tpersp = persp{ii};
        
        % switch statement to define the plots perspective
        switch tpersp
            case{'lh-sag'}
                viewCoords = [-90, 0];
                indx = [slice(ii) 0 0];
                xbnd = [ xbnd; indx ];
            case{'rh-sag'}
                viewCoords = [90, 0];
                indx = [slice(ii) 0 0];
                xbnd = [ xbnd; indx ];
            case{'axial'}
                viewCoords = [0, 90];
                indx = [0 0 slice(ii)];
                xbnd = [ xbnd; indx ];
            case{'coronal'}
                viewCoords = [180, 0];
                indx = [0 slice(ii) 0];
                xbnd = [ xbnd; indx ];
            otherwise
                error('Invalid slice view requested.');
        end
        
        % plot the brain slice
        fh = mbaDisplayBrainSlice(anat, indx, gca); hold on;
        
    end
    
else
    
    % set the defaults w/o slices
    viewCoords = [ -37.5, 30 ];
    indx = [];
    xbnd = [ 0 0 0 ];
    
end

% switch statement for which connections to model
switch conn
    
    case{'path', 'tract'}
        
        for tct = 1:length(tract)
            [fh, light_h] = mbaDisplayConnectome(fgOut{tct}.fibers, fh, color{tct}, 'single');
            delete(light_h);
        
        end
        
        % if cropping is set
        if crop
            xyz = initView(fgOut, xbnd);
        else
            xyz = nan(1, 6);
        end
        
    case{'neighborhood'}
        
        for tct = 1:length(tract)
            
            [fh, light_h] = mbaDisplayConnectome(pnOut{tct}.fibers, fh, colors{2}, 'single');
            delete(light_h);
        
        end
        
        if crop
            xyz = initView(pnOut, xbnd);
        else
            xyz = nan(1, 6);
        end
        
    case{'both'}
        
        for tct = 1:length(tract)
            
            [fh, light_h1] = mbaDisplayConnectome(fgOut{tct}.fibers, fh, color{tct}, 'single');
            delete(light_h1);
            
            [fh, light_h2] = mbaDisplayConnectome(pnOut{tct}.fibers, fh, colors{2}, 'single');
            delete(light_h2);
            
        end
        
        if crop
            xyz = initView(pnOut, xbnd);
        else
            xyz = nan(1, 6);
        end
    
end

% plot this way every time
if crop
    axis(xyz);
end

view(viewCoords);
lh = camlight('right');
lighting phong

end

% create axis call that shows all streamlines / slices requested
function [ out ] = initView(fibers, xbnd)

% find min / max requested x/y/z boundaries
sbnd = minmax(xbnd');

% combine all the nodes of all the streamlines
fibs = [];
for ii = 1:length(fibers)
    ftct = fibers{ii}.fibers;
    for jj = 1:length(ftct)
        fibs = [ fibs ftct{jj} ];
    end
end

%  find min / max [x,y,z]
axCoords = minmax(fibs);

%  compute 10% padding of [x,y,z]
axPadding = abs(axCoords(:,1) - axCoords(:,2)) * 0.10;

% apply default limits to axis
axOut = [ round(axCoords(:,1) - axPadding), round(axCoords(:,2) + axPadding) ];

% adjust limits based on requested slices so they appear
for kk = 1:3
    
    % pull slice indices
    lims = sbnd(kk,:);
    fbnd = axOut(kk,:);
    
    % if no requested slices, skip the check
    if all(lims == 0)
        continue
    end
    
    % pull min and max values and indices
    % this is how you generalize left and right
    [ mn, ni ] = min(fbnd);
    [ mx, xi ] = max(fbnd);
    
    % if the minimum voxel is lower than the minimum slice
    if mn > min(lims)
        % lower the minimum
        fbnd(ni) = min(lims);
    end
    
    % if the maximum voxel is higher than the maximum slice
    if mx < max(lims)
        % raise the maximum
        fbnd(xi) = max(lims);
    end
    
    % store output
    axOut(kk,:) = fbnd;
    
end

% create vector to pass to axis
out = reshape(axOut', 1, 6);

end



