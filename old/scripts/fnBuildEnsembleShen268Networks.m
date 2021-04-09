function [ pconn, rois, omat, olab, files ] = fnBuildEnsembleShen268Networks(subj, ens, nclust, cacheDir)
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

% built less awful names for fe structures - no conditional names
files.fe  = ['fe_' subj '_ensemble_' ens '_REP01_500000.mat'];
files.out = ['fn_' subj '_Shen268_ensemble_' ens '_rep01.mat'];

% build paths files to save
%files.fedir  = '/N/dc2/projects/lifebid/lifeconn';
files.path    = '/N/dc2/projects/lifebid/lifeconn';
files.fedir   = fullfile(files.path, 'subjects', subj, 'fibers');
files.favol   = fullfile(files.path, 'subjects', subj, 'micro', 'fa.nii.gz');
files.outpath = fullfile(files.path, 'subjects', subj, 'networks');
files.output  = fullfile(files.outpath, files.out);
files.labels  = fullfile(files.path, 'subjects', subj, 'anat', [ subj '_shen268_inflated.nii.gz' ]);

%% load files

display('Loading data...');

% load fe structure
load(fullfile(files.fedir, files.fe));

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

% workspace cleaning
save(files.output, 'rois');
clear parc rois

% clean networks
pconn = feCleanPairedConnections(fg, pconn, 'all');
pconn = feCleanPairedConnections(fg, pconn, 'nzw');

% compute tract profiles
favol = niftiRead(files.favol);

pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', favol, 'fa');
pconn = feTractProfilePairedConnections(fg, pconn, 'nzw_clean', favol, 'fa');

clear favol

% virtual lesion matrix
pconn = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, 'nzw');
pconn = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, 'nzw_clean');

%% create adjacency matrices

% create all streamline matrices
[ amat, alab ] = feCreateAdjacencyMatrices(pconn, 'all');
[ nmat, nlab ] = feCreateAdjacencyMatrices(pconn, 'nzw');
[ zmat, zlab ] = feCreateAdjacencyMatrices(pconn, 'nzw_clean');

% combine outputs and labels into 1 matrix
omat = cat(3, amat, nmat, zmat);
olab = [ alab nlab zlab ];

% add stats once I am sure on order
[ glob{1}, node{1}, nets{1} ] = fnNetworkStats(omat(:,:,2), 16);
[ glob{2}, node{2}, nets{2} ] = fnNetworkStats(omat(:,:,10), 16);
[ glob{3}, node{3}, nets{3} ] = fnNetworkStats(omat(:,:,18), 16);

%% save whole output

save(files.output, 'pconn', 'omat', 'olab', 'glob', 'node', 'nets', 'files', '-append');

% remove parallel pool
delete(pool);

end
