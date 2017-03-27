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

% due to alignment frustrations, just plot a binary mask of tract

% create (x, y, z) imagespace of tract
fibNodes = fe.fg.fibers(ind_tract);
fibNodes = cell2mat(fibNodes(:)')';
fibNodes = round(fibNodes) + 1;
fibNodes = unique(fibNodes, 'rows');

% look at feGet
% there may just be a call to that for the indices of what I want...

% % create output image
% imgOut = zeros(imgsize);
% 
% % store wvl and ovl data into volumes
% for ii = 1:size(fibNodes, 1)
%     imgOut(fibNodes(ii, 1), fibNodes(ii, 2), fibNodes(ii, 3)) = 1;
% end

% % create output nifti
% VLmap = niftiCreate('data', imgOut, 'fname', outfile, ...
%                    'qto_xyz', fe.life.xform.img2acpc, ...
%                    'qto_ijk', fe.life.xform.acpc2img);
%                
% % save file
% niftiWrite(VLmap, outfile);

%% plotting VL data

% THESE NEXT 2 LINES ARE HOW TO DO THE VOXEL INTERSECTION FOR LINK NETWORK FROM TENSOR

% pull subtensor of the connection
[ inds, ~ ] = find(fe.life.M.Phi(:, :, ind_tract));

% pull the unique voxels of the connection
voxel_ind = unique(inds(:, 2));

% I can pull precomputed virtual lesion data from pconn
% I just have to convert tensor indices to voxel coordinates

% voxel indices are pulled from converted fiber indices, not the whole volume
% do we even have a mapping back to image space w/ this?
imgsize = fe.life.imagedim(1:3);

% re-create tensor indices b/c they are not the same as image indices
%fibers = fe.fg.fibers;
fibers = fgGet(fe, 'fibers acpc');
fibers = cell2mat(fibers(:)'); 
voxel_coord = ceil(fibers) + 1;
cols = sub2ind(imgsize, voxel_coord(1,:)', voxel_coord(2,:)', voxel_coord(3,:)');
% this is how the voxel indices are originally created

% pull indices from cols
voxCoords = cols(voxel_ind);
% these column indices index back into a rounded fg structure, not the image
% is there even a mapping from these coordinates to image space?

% % convert tensor indices into voxel coordinates
[ voxi, voxj, voxk ] = ind2sub(imgsize, voxCoords);
%[ voxi, voxj, voxk ] = ind2sub(imgsize, voxel_ind); % all coordinates are somewhere along bottom slice
voxCoords = [ voxi voxj voxk ];
clear voxi voxj voxk

% create empty image of zeros that matches coordinate space
wvlimg = zeros(imgsize);
ovlimg = zeros(imgsize);

% store wvl and ovl data into volumes
for ii = 1:size(outCoords, 1)
    wvlimg(voxCoords(ii, 1), voxCoords(ii, 2), voxCoords(ii, 3)) = z.lesion.rmse.all(ii);
    ovlimg(voxCoords(ii, 1), voxCoords(ii, 2), voxCoords(ii, 3)) = z.nolesion.rmse.all(ii);
end

% compute difference between wvl, ovl images
dvlimg = ovlimg - wvlimg;

% merge images into 1 4D image output
imgOut = cat(4, wvlimg, ovlimg, dvlimg);

% create output nifti
VLmap = niftiCreate('data', imgOut, 'fname', outfile, ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);
% any additional parameters?
               
% save file
niftiWrite(VLmap, outfile);

end

