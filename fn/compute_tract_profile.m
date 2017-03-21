function [ fgst, fh, fg ] = compute_tract_profile(fe, pconn, tract, cleaned, msvol, color)
% compute_tract_profile retruns a computed tract profile for precomputed edge
%     using the network outputs, provid a microstructural map in acpc space to compute     
%     properties along the course of a particular connection. 
%     Returns the data itself and a figure handle for simple plot
%
% famp = niftiRead('/N/dc2/projects/lifebid/2t1/HCP/105115/fibers_new/dwi_data_b2000_aligned_trilin_fa.nii.gz');
%
% Brent McPherson (c), 2017
%

% assemble name from pconn structure
roi1 = strrep(pconn{tract}.roi1, '_label', '');
roi1 = strrep(roi1, '_', '.');
roi2 = strrep(pconn{tract}.roi2, '_label', '');
roi2 = strrep(roi2, '_', '.');
fgName = ['Index: ' num2str(tract) '; ' roi1 '-to-' roi2];

% clean tracks if asked to
switch cleaned
    case{'all'}
        fg_tract = fgCreate('name', fgName, 'colorRgb', color, 'fibers', {fe.fg.fibers{pconn{tract}.end.fibers}}');
    case{'nzfibs'}
        fg_tract = fgCreate('name', fgName, 'colorRgb', color, 'fibers', {fe.fg.fibers{pconn{tract}.end.nzfibs}}');
    case{'cln-all'}
        fg_tract = fgCreate('name', fgName, 'colorRgb', color, 'fibers', {fe.fg.fibers{pconn{tract}.cln.fibers}}');
    case{'cln-nzf'}       
        fg_tract = fgCreate('name', fgName, 'colorRgb', color, 'fibers', {fe.fg.fibers{pconn{tract}.cln.nzfibs}}');
end

display(['Computing profile for ' fgName]);

% transform to acpc space
fg = dtiXformFiberCoords(fg_tract, fe.life.xform.img2acpc, 'acpc');

% find stats on FA map w/ vistasoft
fgst = dtiComputeDiffusionPropertiesAlongFG(fg, msvol, [], [], 100);

%% plot - move to a separate fxn later

% compute y axis labels - low / middle / high point
uplw = minmax(fgst');
mdpt = ((max(fgst)-min(fgst)) / 2) + min(fgst);
ybnd = [uplw(1) mdpt uplw(2)];

% plot the FA tract profile
fh = figure('Position', [900 780 875 450]); hold on;
plot(fgst, 'LineWidth', 4, 'Color', color);
title(fgName);
set(gca, 'XTick', [0:50:100], 'YTick', ybnd);
xlabel('Nodes Along Fiber');
ylabel('FA Value');

end

