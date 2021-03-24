%% load development data
%

load('/geode2/home/u010/bcmcpher/Carbonate/fine-data/output_fe.mat');

fg = feGet(fe, 'fibersacpc');
wght = feGet(fe, 'fiberweights');

nTheta = feGet(fe, 'nbvals');
M = feGet(fe, 'model');
dsig = feGet(fe, 'dsigdemeaned by voxel');
S0 = feGet(fe, 'b0signalimage');
%feROI = feGet(fe, 'roicoords');

clear fe

% load full res anatomical
anat = niftiRead('/geode2/home/u010/bcmcpher/Carbonate/fine-data/t1.nii.gz');

% load labels in diffusion space
%fsInflateDK('/geode2/home/u010/bcmcpher/Carbonate/fine-data/aparc+aseg_dwi.nii.gz', 2, 'vert', false, '/geode2/home/bcmcpher/Carbonate/fine-data/parc.nii.gz');
node = niftiRead('/geode2/home/u010/bcmcpher/Carbonate/fine-data/parc.nii.gz');
[ lname, lobe ] = fsMergeDK('/geode2/home/u010/bcmcpher/Carbonate/fine-data/parc.nii.gz', '/geode2/home/u010/bcmcpher/Carbonate/fine-data/lobes.nii.gz');
%parc = niftiRead('/geode2/home/u010/bcmcpher/Carbonate/fine-data/lobes.nii.gz');

fa = niftiRead('/geode2/home/u010/bcmcpher/Carbonate/fine-data/dt6_fa.nii.gz');
dt6 = dtiLoadDt6('/geode2/home/u010/bcmcpher/Carbonate/fine-data/dti/dt6.mat');

nname = {'ctx-lh-bankssts', 'ctx-lh-caudalanteriorcingulate', 'ctx-lh-caudalmiddlefrontal', ...
         'ctx-lh-cuneus', 'ctx-lh-entorhinal', 'ctx-lh-fusiform', 'ctx-lh-inferiorparietal', ...
         'ctx-lh-inferiortemporal', 'ctx-lh-isthmuscingulate', 'ctx-lh-lateraloccipital', ...
         'ctx-lh-lateralorbitofrontal', 'ctx-lh-lingual', 'ctx-lh-medialorbitofrontal', ...
         'ctx-lh-middletemporal', 'ctx-lh-parahippocampal', 'ctx-lh-paracentral', ...
         'ctx-lh-parsopercularis', 'ctx-lh-parsorbitalis', 'ctx-lh-parstriangularis', ...
         'ctx-lh-pericalcarine', 'ctx-lh-postcentral', 'ctx-lh-posteriorcingulate', ...
         'ctx-lh-precentral', 'ctx-lh-precuneus', 'ctx-lh-rostralanteriorcingulate', ...
         'ctx-lh-rostralmiddlefrontal', 'ctx-lh-superiorfrontal', 'ctx-lh-superiorparietal', ...
         'ctx-lh-superiortemporal', 'ctx-lh-supramarginal', 'ctx-lh-frontalpole', ...
         'ctx-lh-temporalpole', 'ctx-lh-transversetemporal', 'ctx-lh-insula', ...
         'ctx-rh-bankssts', 'ctx-rh-caudalanteriorcingulate', 'ctx-rh-caudalmiddlefrontal', ...
         'ctx-rh-cuneus', 'ctx-rh-entorhinal', 'ctx-rh-fusiform', 'ctx-rh-inferiorparietal', ...
         'ctx-rh-inferiortemporal', 'ctx-rh-isthmuscingulate', 'ctx-rh-lateraloccipital', ...
         'ctx-rh-lateralorbitofrontal', 'ctx-rh-lingual', 'ctx-rh-medialorbitofrontal', ...
         'ctx-rh-middletemporal', 'ctx-rh-parahippocampal', 'ctx-rh-paracentral', ...
         'ctx-rh-parsopercularis', 'ctx-rh-parsorbitalis', 'ctx-rh-parstriangularis', ...
         'ctx-rh-pericalcarine', 'ctx-rh-postcentral', 'ctx-rh-posteriorcingulate', ...
         'ctx-rh-precentral', 'ctx-rh-precuneus', 'ctx-rh-rostralanteriorcingulate', ...
         'ctx-rh-rostralmiddlefrontal', 'ctx-rh-superiorfrontal', 'ctx-rh-superiorparietal', ...
         'ctx-rh-superiortemporal', 'ctx-rh-supramarginal', 'ctx-rh-frontalpole', ...
         'ctx-rh-temporalpole', 'ctx-rh-transversetemporal', 'ctx-rh-insula'}';

%length = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), fg.fibers, 'UniformOutput', true);

clc
