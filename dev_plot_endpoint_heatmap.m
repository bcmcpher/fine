%% end point density plot

load('test/fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat');

% reload image space ep - save those down
% initialize endpoint outputs
iep1 = zeros(length(fe.fg.fibers), 3);
iep2 = zeros(length(fe.fg.fibers), 3);

% for every fiber, pull the end points
for ii = 1:length(fe.fg.fibers)
    iep1(ii,:) = fe.fg.fibers{ii}(:,1)';
    iep2(ii,:) = fe.fg.fibers{ii}(:,end)';
end

% round all the endpoints
iep1 = round(iep1) + 1;
iep2 = round(iep2) + 1;

% combine and reconvert fiber endpoints to image space
ep = [ iep1; iep2 ];

% sanity check of low/high dims
minmax(ep')

% create empty data space
img = zeros(fe.life.imagedim(1:3));
% THIS IS NOT THE SAME SIZE AS APARC
% IN REGISTER, BUT NOT THE SAME SHAPE

% does this represent the translation offset?
aparc.qto_xyz - fe.life.xform.img2acpc;

% for every endpoint, create a density map of rois
for ii = 1:length(ep)
    img(ep(ii,1), ep(ii, 2), ep(ii, 3)) = img(ep(ii,1), ep(ii, 2), ep(ii, 3)) + 1;
end

% create output .nii
hmap = niftiCreate('data', img, 'fname', 'hmap/fe_out.nii.gz', ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);
               
% save file
niftiWrite(hmap, hmap.fname);


