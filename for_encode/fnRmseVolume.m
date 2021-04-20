function [ rmse ] = fnRmseVolume(fe, smooth, reslice, outfile, norm)
%[ rmse ] = fnRmseVolume(fe, smooth, outfile, norm);
%    This function creates a RMSE volume from a fit fe structure.
%    Will optionally reslice to a specified volume and / or smooth. 
%    It returns a nifti object with the extracted data.
% 
% INPUTS:
%     fe      - an evaluated fiber evaluation structure (LiFE has been estimated)
%     smooth  - perform 3D smoothing of provided kernel size 
%     reslice - the nifti file the RMSE will be resliced to match
%     outfile - (optional) the name of the .nii.gz file to write the volume to
%     norm    - correct (normalize) the RMSE based on the b0 (default = true)
% 
% OUTPUTS:
%     rmse    - a .nii structure containing the computed heatmap data
%
% Brent McPherson (c) 2017 Indiana University
%

% use normalized RMSE by default
if(~exist('norm', 'var') || isempty(norm))
    norm = 1;
end

% turn off smoothing if it's not requested
if(~exist('smooth', 'var') || isempty(smooth))
   smooth = [];
end

% deal with reslicing input
if(~exist('reslice', 'var') || isempty(reslice)) % if it doesn't exist, skip
    reslice = [];
else
    if ischar(reslice) && ~isfield(reslice, 'nifti_type') % try to load if it's a string and not a nifti already
        try
            reslice = niftiRead(reslice);
        catch
            warning('Unable to load nifti image for reslicing output. Defaulting back to fe space for output.');
            reslice = [];
        end
    else % warn if it can't load the input as a nifti or if the input does not have the nifti_type field
        warning('Unknown input passed for reslice option. Defaulting back to fe space for output.');
    end
end

% set outfile to empty if a filename isn't provided (will not write)
if(~exist('outfile', 'var') || isempty(outfile))
   outfile = [];
end

% create output dimensions
out_xyz = fe.life.xform.img2acpc;
out_ijk = fe.life.xform.acpc2img;
out_dim = fe.life.imagedim(1:3);

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

% create empty data space
img = zeros(out_dim);

% for every roi coordinate pull the rmse
for ii = 1:size(roi, 1)
    img(roi(ii, 1), roi(ii, 2), roi(ii, 3)) = vals(ii);
end

% create output nifti in fe space
rmse = niftiCreate('data', img, 'fname', outfile, ...
                   'qto_xyz', out_xyz, ...
                   'qto_ijk', out_ijk);
% any additional parameters?

% if reslice is not empty, reslice the rmse nifti to requested space
if ~isempty(reslice)
    rmse = mrAnatResampleToNifti(rmse, reslice);
end

% smooth nifti data if requested
if ~isempty(smooth)
    disp(['Smoothing output with a [ ' num2str(smooth) ' ] gaussian kernel...']);
    rmse.data = smooth3(rmse.data, 'gaussian', smooth);
end

% save file if outfile is defined
if ~isempty(outfile)
    niftiWrite(rmse, outfile);
end

end
