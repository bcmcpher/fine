function [ VLmap ] = fnPlotVirtualLesion(fe, pconn, tract, label, outfile)
%fnPlotVirtualLesion() creates a volume of a connections change in signal
% estimated by feComputeVirtualLesion
%   
% Brent McPherson, (c) 2017
%

% pull streamline indices for a particular connection
ind_tract = getfield(pconn{tract}, label, 'indices');

% dummy check if the connection is empty
if isempty(ind_tract)
    error('This connection has no streamlines. There is no Virtual Lesion output.');
end

%% plotting VL data

% pull subtensor of the connection
[ inds, ~ ] = find(fe.life.M.Phi(:, :, ind_tract));

% pull the unique voxels of the connection
voxel_ind = unique(inds(:, 2));

% pull image coordinates of tensor indexed voxels
coords = feGet(fe, 'coordsfromfibers', ind_tract);

% pull the voxelwise rmse coords for the tract
vals = feGet(fe, 'voxrmse', coords);

% make sure this is the right value...

% total R^2 explained
%z = feGet(fe, 'totalr2');
%z = feGet(fe, 'totpve');

% broke...
%z = feGet(fe, 'voxelvarianceexplainedvoxelwise', fibNodes); % i think this will return the whole volume

% create image size
imgsize = fe.life.imagedim(1:3);

% create empty image of zeros that matches coordinate space
imgOut = zeros(imgsize);

% store wvl and ovl data into volumes
for ii = 1:size(coords, 1)
    imgOut(coords(ii, 1), coords(ii, 2), coords(ii, 3)) = vals(ii);
end

% create output nifti
VLmap = niftiCreate('data', imgOut, 'fname', outfile, ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);
% any additional parameters?
               
% save file
niftiWrite(VLmap, outfile);

end

