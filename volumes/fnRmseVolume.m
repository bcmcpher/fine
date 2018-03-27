function [ rmse ] = fnRmseVolume(fe, outfile, norm, smooth)
%rmse = fnRmseVolume() creates a RMSE volume from an fe structure.
% Will optionally reslice to a specified volume and / or smooth. 
% It returns a nifti object with the extracted data.
% 
% INPUTS:
%     fe      - an evaluated fiber evaluation structure (LiFE has been performed)
%     outfile - the name of the .nii.gz file to write the volume to
%     norm    - whether or not to correct RMSE based on b0
%     smooth  - perform 3D smoothing of provided kernel size
% 
% OUTPUTS:
%     rmse    - a .nii structure containing the computed heatmap data
%
% TODO: 
% - figure out reslice option
%
% Brent McPherson (c) 2017 Indiana University
%


% use normalized RMSE by default
if(~exist('norm', 'var') || isempty(norm))
    norm = 1;
end

% % if no new space is provided, use the fe structure for output size
% if(~exist('reslice', 'var') || isempty(reslice))
%     disp('Image will be sliced in fe space...');
%     do_reslice = 0;
%     out_xyz = fe.life.xform.img2acpc;
%     out_ijk = fe.life.xform.acpc2img;
%     out_dim = fe.life.imagedim(1:3);
% 
% % otherwise use the data from the new space for output size
% else
%     disp(['Image will be resliced into ' reslice.fname ' space...']);
%     do_reslice = 1;
%     out_xyz = reslice.qto_xyz;
%     out_ijk = reslice.qto_ijk;
%     out_dim = size(reslice.data);    
% end

% create output dimensions
out_xyz = fe.life.xform.img2acpc;
out_ijk = fe.life.xform.acpc2img;
out_dim = fe.life.imagedim(1:3);

% turn off smoothing if it's not requested
if(~exist('smooth', 'var') || isempty(smooth))
   smooth = [];
end

% % default smoothing if data is resliced and no smoothing is requested 
% if (do_reslice == 1 && isempty(smooth))
%     smooth = [ 3 3 3 ];
% end

% decide how to compute RMSE and get it
switch norm
    case {1, 'norm'}
        disp('Computing normalized RMSE...');
        vals = feGet(fe, 'voxrmses0norm');
    otherwise
        disp('Computing RMSE...');
        vals = feGet(fe, 'voxrmse');
end

% extract VOI
roi = feGet(fe, 'roicoords');

% % convert VOI to acpc space
% roi_ac = mrAnatXformCoords(fe.life.xform.img2acpc, roi);
% 
% % move VOI back into output space
% roi = round(mrAnatXformCoords(out_ijk, roi_ac));

% create empty data space
img = zeros(out_dim);

% for every endpoint, create a density map of rois
for ii = 1:size(roi, 1)
    img(roi(ii, 1), roi(ii, 2), roi(ii, 3)) = vals(ii);
end

% smooth data if requested
if ~isempty(smooth)
    disp(['Smoothing output with a [ ' num2str(smooth) ' ] gaussian kernel...']);
    img = smooth3(img, 'gaussian', smooth);
end

% replace NaNs?

% create output nifti
rmse = niftiCreate('data', img, 'fname', outfile, ...
                   'qto_xyz', out_xyz, ...
                   'qto_ijk', out_ijk);
% any additional parameters?
               
% save file
niftiWrite(rmse, outfile);

end

