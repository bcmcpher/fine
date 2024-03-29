function [ pconn, rois, omat, olab, files ] = fnBuildShen268Networks(subj, track, lmax, rep, nclust, cacheDir)
% fnBuildShen268Networks creates all the data I need for the paper from a fe structure.
%   This will probably get renamed because this is a short, simple process
%   to run that is canonical with many pipelines.
%
% INPUTS: all inputs are strings
%
%   subj  - subject ID, matches folder name
%   track - tracking parameter, 'SD_STREAM', 'SD_PROB', 'tensor_'
%   lmax  - lmax parameter
%   rep   - repetition number
%
% Example usage:
%   
% subj = '105115'; track = 'SD_PROB'; lmax='10'; rep='01'; nclust = 12; cacheDir='/N/dc2/projects/lifebid/lifeconn/matlab/cache';
% [ pconn, rois, omat, olab, files ] = fnBuildShen268Networks('105115', 'SD_PROB', '10', '01', 12, '/N/dc2/projects/lifebid/lifeconn/matlab/cache');
%

%% build paths

% build file names if tensor is called
if strcmp(track, {'tensor_'})
    files.fe = ['fe_structure_' subj '_STC_run01_' track '_lmax' lmax '_connNUM' rep '.mat'];
    files.out     = ['fn_' subj '_Shen268_' track '_rep' rep '.mat'];
    files.csv     = ['fn_' subj '_Shen268_' track '_rep' rep ];
else
    files.fe = ['fe_structure_' subj '_STC_run01_' track '_lmax' lmax '_connNUM' rep '.mat'];
    files.out     = ['fn_' subj '_Shen268_' track '_lmax' lmax '_rep' rep '.mat'];
    files.csv     = ['fn_' subj '_Shen268_' track '_lmax' lmax '_rep' rep ];
end

% build paths files to save
files.fedir  = '/N/dc2/projects/lifebid/code/ccaiafa/Caiafa_Pestilli_paper2015/Revision_Feb2017/Results/Single_TC/fe_structures';
files.path    = '/N/dc2/projects/lifebid/lifeconn';
files.favol   = fullfile(files.path, 'subjects', subj, 'micro', 'fa.nii.gz');
files.outpath = fullfile(files.path, 'subjects', subj, 'networks');
files.output  = fullfile(files.outpath, files.out);
files.labels  = fullfile(files.path, 'subjects', subj, 'anat', [ subj '_shen268_inflated.nii.gz' ]);

% conditionals to build subject names, because that's not annoying

% data sets with funny names somewhere in the process

% because these FE file names are still different by data set...

STN = {'FP', 'KK', 'HT', 'MP'};
HC7 = {'7t_111312', '7t_125525', '7t_102311', '7t_109123'};

% stanford data has this extra file info after subject ID in the fe structure names, not paths
if any(strcmp(subj, STN))
    nsubj = [ subj '_96dirs_b2000_1p5iso' ];
    files.fe = ['fe_structure_' nsubj '_STC_run01_' track '_lmax' lmax '_connNUM' rep '.mat'];
end

% 7T subjects are 7t_*, except in the fe structure names, where they're just the number
if any(strcmp(subj, HC7))
    nsubj = strrep(subj, '7t_', '');
    files.fe = ['fe_structure_' nsubj '_STC_run01_' track '_lmax' lmax '_connNUM' rep '.mat'];
end

clear STN HC7

%% load files

display('Loading data...');

% load fe structure
load(fullfile(files.fedir, subj, files.fe));

% extract all needed out of FE
fg               = feGet(fe,   'fg acpc'); 
fascicle_length  = fefgGet(fg, 'length');
fascicle_weights = feGet(fe,   'fiber weights');
nTheta           = feGet(fe,   'nbvals');
M                = feGet(fe,   'model');
measured_dsig    = feGet(fe,   'dsigdemeaned by voxel');
clear fe

% load labeled aparc+aseg volume
parc = niftiRead(files.labels);

%% start parallel pool

display(['Opening parallel pool with ', num2str(nclust), ' cores...']);

% create parallel cluster object
clust = parcluster;

% set number of cores from arguments
clust.NumWorkers = nclust;

% set temporary cache directory
tmpdir = tempname(cacheDir);

% make cache dir
OK = mkdir(tmpdir);

% check and set cachedir location
if OK
    % set local storage for parpool
    clust.JobStorageLocation = tmpdir;
end

% start parpool - close parpool at end of fxn
pool = parpool(clust, nclust, 'IdleTimeout', 120);

clear tmpdir OK

% clear subj track lmax rep nclust cacheDir

%% create networks

% assign streamline endpoints to labeled volume

% Run operations
[ pconn, rois ] = feCreatePairedConnections(parc, fg.fibers, fascicle_length, fascicle_weights);

% clean networks
%pconn = feCleanPairedConnections(fg, pconn, 'all');
%pconn = feCleanPairedConnections(fg, pconn, 'nzw');

% compute tract profiles
favol = niftiRead(files.favol);

pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', favol, 'fa');
%pconn = feTractProfilePairedConnections(fg, pconn, 'nzw_clean', favol, 'fa');

clear favol

% virtual lesion matrix
pconn = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, 'nzw');
%pconn = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, 'nzw_clean');

%% create adjacency matrices

% create all streamline matrices
[ amat, alab ] = feCreateAdjacencyMatrices(pconn, 'all');
[ nmat, nlab ] = feCreateAdjacencyMatrices(pconn, 'nzw');
%[ zmat, zlab ] = feCreateAdjacencyMatrices(pconn, 'nzw_clean');

% combine outputs and labels into 1 matrix
omat = cat(3, amat, nmat);
olab = [ alab nlab ];

% add stats once I am sure on order
[ glob{1}, node{1}, nets{1} ] = fnNetworkStats(omat(:,:,2), 16);
[ glob{2}, node{2}, nets{2} ] = fnNetworkStats(omat(:,:,10), 16);
%[ glob{3}, node{3}, nets{3} ] = fnNetworkStats(omat(:,:,18), 16);

%% save whole output

save(files.output, 'pconn', 'rois', 'omat', 'olab', 'glob', 'node', 'nets', 'files', '-v7.3');

% remove parallel pool
delete(pool);

% %% save what O3D needs
% 
% % create labels vector
% labels = zeros(length(rois), 1);
% for ii = 1:length(rois)
%     labels(ii) = rois{ii}.label;
% end
% 
% % write .csv files
% dlmwrite(fullfile(files.outpath, [ files.csv '_count.csv' ]), omat(:,:,1));
% dlmwrite(fullfile(files.outpath, [ files.csv '_density.csv' ]), omat(:,:,2));
% dlmwrite(fullfile(files.outpath, [ files.csv '_length.csv' ]), omat(:,:,3));
% dlmwrite(fullfile(files.outpath, [ files.csv '_labels.csv' ]), labels);

end
