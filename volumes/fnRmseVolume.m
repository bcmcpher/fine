function [ rmse ] = fnRmseVolume(fe, outfile, norm)
%% compute rmse volume
%
% Create RMSE map - matches ACPC alignment, in register w/ anat, in DWI space
%
% TODO:
% - other normalization options?
% - reslice options?
% - make sure indices line up
%

% use regular RMSE by default
if(~exist('norm', 'var') || isempty(norm))
    norm = 0;
end

% decide how to compute RMSE and get it
switch norm
    case {1, 'norm'}
        display('Computing normalized RMSE...');
        vals = feGet(fe, 'voxrmses0norm');
    otherwise
        display('Computing RMSE...');
        vals = feGet(fe, 'voxrmse');
end

% pull ROI coordinates
roi = feGet(fe, 'roicoords');

% create empty data space
img = zeros(fe.life.imagedim(1:3));

% for every endpoint, create a density map of rois
for ii = 1:size(roi, 1)
    img(roi(ii, 1), roi(ii, 2), roi(ii, 3)) = vals(ii);
end

% create output nifti
rmse = niftiCreate('data', img, 'fname', outfile, ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);
% any additional parameters?
               
% save file
niftiWrite(rmse, outfile);

end

