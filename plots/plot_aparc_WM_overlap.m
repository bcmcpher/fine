%% quick plot of ROI voxels and voxels used as endpoint white matter intersections
% Brent McPherson
% 20170116
%

%% load data

% load an fe structure for a subject
load fe_structure_105115_STC_run01_SD_PROB_lmax10_connNUM01.mat

aparc = niftiRead('aparc+aseg.nii.gz');
aparc.data(aparc.data < 1000) = 0;
aparc.data(aparc.data == 1000 | aparc.data == 2000) = 0;

iparc = niftiRead('inflated_labels.nii.gz');

% % load the individual freesurfer ROI to check the end points
% roi1 = dtiImportRoiFromNifti('rois/rh_cuneus_label.nii.gz');
% roi2 = dtiImportRoiFromNifti('rois/rh_lateraloccipital_label.nii.gz');
% roi3 = dtiImportRoiFromNifti('rois/rh_lingual_label.nii.gz');

%% apply transform to points, find indices within ROI of voxels

% % move ROIs to acpc space
% vx{1}.coords = mrAnatXformCoords(fe.life.xform.acpc2img, roi1.coords);
% vx{2}.coords = mrAnatXformCoords(fe.life.xform.acpc2img, roi2.coords);
% vx{3}.coords = mrAnatXformCoords(fe.life.xform.acpc2img, roi3.coords);

% create roi of aparc data
[x, y, z ] = ind2sub(size(aparc.data), find(aparc.data));
aparc_mask = [ x, y, z ];
clear x y z

[x, y, z ] = ind2sub(size(iparc.data), find(iparc.data));
iparc_mask = [ x, y, z ];
clear x y z

% move image coordinates to acpc coordinates
% also try xform: aparc.qto_xyz
vxa = mrAnatXformCoords(fe.life.xform.img2acpc, aparc_mask);
roi = mrAnatXformCoords(fe.life.xform.img2acpc, fe.roi.coords);
ixa = mrAnatXformCoords(fe.life.xform.img2acpc, iparc_mask);

% % fix round off
% vx{1}.coords = ceil(vx{1}.coords);
% vx{2}.coords = ceil(vx{2}.coords);
% vx{3}.coords = ceil(vx{3}.coords);

%vxar = ceil(vxa);
%vxar = floor(vxa + 0.50);
vxar = round(vxa) + 1;

% % find intersect of cortical ROI in LiFE white matter ROI
% index{1} = ismember(fe.roi.coords, vx{1}.coords, 'rows');
% index{2} = ismember(fe.roi.coords, vx{2}.coords, 'rows');
% index{3} = ismember(fe.roi.coords, vx{3}.coords, 'rows');

indx = ismember(fe.roi.coords, vxar, 'rows');

% % proportion of ROI used as end points
% sum(index{1}) / size(vx{1}.coords, 1)
% sum(index{2}) / size(vx{2}.coords, 1)
% sum(index{3}) / size(vx{3}.coords, 1)

sum(indx) / size(vxar, 1)

% % find coordinates of ROI in fe.roi
% feROIint{1} = intersect(fe.roi.coords, vx{1}.coords, 'rows');
% feROIint{2} = intersect(fe.roi.coords, vx{2}.coords, 'rows');
% feROIint{3} = intersect(fe.roi.coords, vx{3}.coords, 'rows');

feROI = intersect(fe.roi.coords, vxar, 'rows');

%% plot the data

% red colors are scaled to line up with the ROI

figure;
hold on;

% plot the 3 rois with their intersection into white matter
% loop instead?
% better color options?
% more optimal fills / shapes to use?

plot3(vxar(:,1), vxar(:,2), vxar(:,3), '.', 'markersize', 10, 'color', 'blue');
plot3(feROI(:, 1), feROI(:, 2), feROI(:, 3), 'o', 'color', [1 0 0], 'markerfacecolor',[1 0 0], 'linewidth', 3);

% plot3(vx{1}.coords(:,1), vx{1}.coords(:,2), vx{1}.coords(:,3), '.', 'markersize', 10, 'color', 'blue');
% plot3(feROIint{1}(:, 1), feROIint{1}(:, 2), feROIint{1}(:, 3), 'o', 'color', [1 0 0], 'markerfacecolor',[1 0 0], 'linewidth', 3);
% 
% plot3(vx{2}.coords(:,1), vx{2}.coords(:,2), vx{2}.coords(:,3), '.', 'markersize', 10, 'color', 'green');
% plot3(feROIint{2}(:, 1), feROIint{2}(:, 2), feROIint{2}(:, 3), 'o', 'color', [0.75 0 0], 'markerfacecolor', [0.75 0 0], 'linewidth', 3);
% 
% plot3(vx{3}.coords(:,1), vx{3}.coords(:,2), vx{3}.coords(:,3), '.', 'markersize', 10, 'color', 'cyan');
% plot3(feROIint{3}(:, 1), feROIint{3}(:, 2), feROIint{3}(:, 3), 'o', 'color', [0.50 0 0], 'markerfacecolor',[0.50 0 0], 'linewidth', 3);

% scale axes and set to coronal view
axis equal; axis square
view(0, 0)


