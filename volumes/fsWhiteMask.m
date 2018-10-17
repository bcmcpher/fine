function [ wm ] = fsWhiteMask(aparc, outfile)
%fsWhiteMask takes an in register aparc image and saves out the white
% matter tracking volume
%   Detailed explanation goes here

% define output file
%wmMaskFile = fullfile(outdir, subjects{isbj}, 'anat', 'wm_mask.nii.gz');

% import and convert freesurfer file
%fs_wm = fullfile(fsdir, subjects{isbj}, 'mri', 'aseg.mgz');
%eval(sprintf('!mri_convert  --out_orientation RAS %s %s', fs_wm, outfile));
wm = niftiRead(aparc);

% define labels to convert into WM mask
invals  = [2 41 16 17 28 60 51 53 12 52 13 18 54 50 11 251 252 253 254 255 10 49 46 7];
%invals  = [251 252 253 254 255];

origvals = unique(wm.data(:));

% convert multiple labels to binary image
fprintf('\n[%s] Converting voxels... ',mfilename);

wmCounter=0; noWMCounter=0;
for ii = 1:length(origvals);
    if any(origvals(ii) == invals)
        wm.data( wm.data == origvals(ii) ) = 1;
        wmCounter=wmCounter+1;
    else
        wm.data( wm.data == origvals(ii) ) = 0;
        noWMCounter = noWMCounter + 1;
    end
end

fprintf('converted %i regions to White-matter (%i regions left outside of WM)\n\n',wmCounter,noWMCounter);

% change file header to outfile
wm.fname = outfile;

% save new nifti
niftiWrite(wm, outfile);

end

