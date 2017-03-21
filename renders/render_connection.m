function [ fh, fg ] = render_connection(fe, pconn, tract, persp, cleaned, conn, pnprc)
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

% path to save down figures
%savedir = './figs/';

%% other arguments...

% add default if pnprc is empty / undefined...

% set minimum length of 10mm
threshold_length = 10;

% generate the blue / red colors for the tracts
colors = {[.1 .25 .65], [.75 .25 .1]};

%% run the analysis

% pull name from pconn based on index
roi1 = strrep(pconn{tract}.roi1, '_label', '');
roi1 = strrep(roi1, '_', '.');
roi2 = strrep(pconn{tract}.roi2, '_label', '');
roi2 = strrep(roi2, '_', '.');
fgName = ['Index: ' num2str(tract) '; ' roi1 '-to-' roi2];

% pick tracts to plot
switch cleaned
    case{'all'}
        fg_tract = fgCreate('name', fgName, 'colorRgb', colors{1}, 'fibers', {fe.fg.fibers{pconn{tract}.end.fibers}}');
        ind1 = pconn{tract}.end.fibers;
    case{'nzfibs'}
        fg_tract = fgCreate('name', fgName, 'colorRgb', colors{1}, 'fibers', {fe.fg.fibers{pconn{tract}.end.nzfibs}}');
        ind1 = pconn{tract}.end.nzfibs;
    case{'cln-all'}
        fg_tract = fgCreate('name', fgName, 'colorRgb', colors{1}, 'fibers', {fe.fg.fibers{pconn{tract}.cln.fibers}}');
        ind1 = pconn{tract}.cln.fibers;
    case{'cln-nzf'}       
        fg_tract = fgCreate('name', fgName, 'colorRgb', colors{1}, 'fibers', {fe.fg.fibers{pconn{tract}.cln.nzfibs}}');
        ind1 = pconn{tract}.cln.nzfibs;
    otherwise
        error(['Invalid type of fibers selected: ' cleaned ' is not an option. Please check you input for the ''clean'' argument.']);
end

display(['Rednering ' fgName]);

% if an emtpy fg is requested
if isempty(ind1)
    error('The connection requested contains no streamlines at this parameter.');
end

% pull path neighborhood indices
ind2 = feGet(fe, 'Path Neighborhood', ind1);

% create fiber group
fg = feGet(fe, 'fibers img');

% diff order of extract / align
xform = feGet(fe, 'img2acpcxform');
fg_xform = dtiXformFiberCoords(fg, xform, 'acpc');

% subset path neighborhood from whole brain fg, clear whole brain fg
fg_pathn = fgExtract(fg_xform, ind2, 'keep');

clear fg

% transform edge coords to acpc space, add to catch in output object
fg{1} = dtiXformFiberCoords(fg_tract, xform, 'acpc');

% drop unused field and create final path neighborhood object before filtering
fg_pathn = rmfield(fg_pathn, 'coordspace');
fg_pnplot = fg_pathn;

% SHOULD BE DONE AS PART OF CLEANING IN CONNECTOMES, NOT PLOTTING
% filter fiber lenghts in path neighborhood
c = 1;
fibers = cell(1, length(fg_pathn.fibers));
for ii = 1:length(fg_pathn.fibers)
    if length(fg_pathn.fibers{ii}) > threshold_length
        fibers{c} = fg_pathn.fibers{ii};
        c = c + 1;
    end
end

fg_pathn.fibers = fibers; 
clear fibers

% pull a random subset of the path neigborhood to minimize the rendering
fibs_indx = randsample(1:length(fg_pathn.fibers), round(length(fg_pathn.fibers) * pnprc));
fg_pnplot.fibers = fg_pnplot.fibers(fibs_indx);

% assign path neighborhood as second 
fg{2} = fg_pnplot;

%% plot the requested connection / view / fg

% switch statement to define the plots perspective
switch persp
    case{'lh-sag'}
        viewCoords = [-90, 0];
        slice = [-1 0 0];
    case{'rh-sag'}
        viewCoords = [90, 0];
        slice = [-1 0 0];
    case{'axial'}
        viewCoords = [0, 90];
        slice = [0 0 -1];
    case{'coronal'}
        viewCoords = [180, 0];
        slice = [0 -1 0];
end

% plot the brain slice
fh = mbaDisplayBrainSlice(anatomy, slice, gca); hold on;

% switch statement for which connections to model
switch conn
    case{'path'}
        [fh, light_h] = mbaDisplayConnectome(fg{1}.fibers, fh, colors{1}, 'single');
        delete(light_h);
    case{'neighborhood'}
        [fh, light_h] = mbaDisplayConnectome(fg{2}.fibers, fh, colors{2}, 'single');
        delete(light_h);
    case{'both'}
        [fh, light_h1] = mbaDisplayConnectome(fg{1}.fibers, fh, colors{1}, 'single');
        delete(light_h1);
        [fh, light_h2] = mbaDisplayConnectome(fg{2}.fibers, fh, colors{2}, 'single');
        delete(light_h2);
end

% plot this way every time
view(viewCoords);
camlight right
lighting phong

% save images? - don't run, this breaks the session
%feSavefig(fhNum, 'verbose', 'yes', 'figName', fig.names{iview}, 'figDir', fullFigureOutDir, 'figType', 'jpg');
%print(fh, '-cmyk', '-painters', '-depsc2', '-tiff', '-r500' , '-noui', './renders/test.eps');

end
