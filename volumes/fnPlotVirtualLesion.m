function [ VLmap ] = fnPlotVirtualLesion(fe, pconn, tract, label, outfile)
%fnPlotVirtualLesion() creates a volume of a connections change in signal
% estimated by feComputeVirtualLesion
%   
% Brent McPherson, (c) 2017
%

%% extract data

% pull streamline indices for a particular connection
indices = getfield(pconn{tract}, label, 'indices');

% dummy check if the connection is empty
if isempty(indices)
    error('This connection has no streamlines. Nothing can be done.');
end

% pull precomputed virtual lesion data
vl = getfield(pconn{tract}, label, 'vl');

% dummy check if virtual lesion is there
if isempty(vl)
    error('No virtual lesion field found. Make sure virtual lesion has been run.');
end

%% create plot of VL data

% pull image coordinates of tensor indexed voxels
coords = feGet(fe, 'coordsfromfibers', indices);

% pull the voxelwise rmse coords for the tract
rmse = feGet(fe, 'voxrmse', coords);

% pull all rmse value for lesioned and unlesioned data
vlo = vl.lesion.rmse.all;
vlw = vl.nolesion.rmse.all;

% find the difference - is this sensical?
vld = vlo - vlw;

% make sure these are the right values...

% create image size
imgsize = fe.life.imagedim(1:3);

% create empty image of zeros that matches coordinate space
imgOut = zeros([ imgsize 4 ]);

% store wvl and ovl data into volumes
for ii = 1:size(coords, 1)
    imgOut(coords(ii, 1), coords(ii, 2), coords(ii, 3), 1) = rmse(ii);
    imgOut(coords(ii, 1), coords(ii, 2), coords(ii, 3), 2) = vlo(ii);
    imgOut(coords(ii, 1), coords(ii, 2), coords(ii, 3), 3) = vlw(ii);
    imgOut(coords(ii, 1), coords(ii, 2), coords(ii, 3), 4) = vld(ii);
end

% create output nifti
VLmap = niftiCreate('data', imgOut, 'fname', outfile, ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);
% any additional parameters?
               
% save file
niftiWrite(VLmap, outfile);

end

%% other useful measures to pull?

% relative estimation error for a particular edge
%x = feGet(fe, 'relativeerror', coords);

% would be virtual lesion if it worked - 'psigfibertest' isn't an feGet option
%x = feGet(fe, 'voxrmsevoxelwise', coords);
%y = feGet(fe, 'voxrmsetest', ind_tract, coords);

% path neighborhood indices
%x = feGet(fe, 'pathneighborhood', coords);

% fiber density of each coordinate - fg needs to be path neighborhood?
%x = fefgGet(fe.fg, 'fiberdensity', coords);

% cell array of unique fiber indices in each voxel (path neighborhood)
%x = fefgGet(fe.fg, 'uniquefibersinvox', coords);

% I should understand Westin Shapes - useful extra parameter?
%dtiComputeWestinShapes

% total R^2 explained
%z = feGet(fe, 'totalr2');
%z = feGet(fe, 'totpve');
%z = feGet(fe, 'totalrmse');

% broke...
%z = feGet(fe, 'voxelvarianceexplainedvoxelwise', fibNodes); % i think this will return the whole volume
