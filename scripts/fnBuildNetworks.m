function [ pconn, rois, omat, olab ] = fnBuildNetworks(fe, labels, fa, output, nclust, cacheDir)
%fnBuildShen268Networks creates all the data I need for the paper from a fe structure.
%   This will probably get renamed because this is a short, simple process
%   to run that is canonical with many pipelines.
%
% INPUTS: all inputs are strings
% subj  - subject ID, matches folder name
% track - tracking parameter, 'SD_STREAM', 'SD_PROB', 'tensor_'
% lmax  - lmax parameter
% rep   - repetition number
%
% subj = '105115'; track = 'SD_PROB'; lmax='10'; rep='01'; nclust = 12; cacheDir='/N/dc2/projects/lifebid/lifeconn/matlab/cache';
% [ pconn, rois, omat, olab, files ] = fnBuildShen268Networks('105115', 'SD_PROB', '10', '01', 12, '/N/dc2/projects/lifebid/lifeconn/matlab/cache');
%

%% load files

display('Loading data...');

% load fe structure
load(fe);

% extract all needed out of FE
fg               = feGet(fe,   'fg acpc'); 
fascicle_length  = fefgGet(fg, 'length');
fascicle_weights = feGet(fe,   'fiber weights');
nTheta           = feGet(fe,   'nbvals');
M                = feGet(fe,   'model');
measured_dsig    = feGet(fe,   'dsigdemeaned by voxel');

clear fe

% load labels volume and FA image
parc  = niftiRead(labels);
favol = niftiRead(fa);

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

%% create networks

% assign streamline endpoints to labeled volume
[ pconn, rois ] = feCreatePairedConnections(parc, fg.fibers, fascicle_length, fascicle_weights);

% clean networks
%pconn = feCleanPairedConnections(fg, pconn, 'all');
%pconn = feCleanPairedConnections(fg, pconn, 'nzw');

% run virtual lesion
pconn = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, 'nzw');

% compute tract profiles
pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', favol, 'fa');
%pconn = feTractProfilePairedConnections(fg, pconn, 'nzw_clean', favol, 'fa');

%% create adjacency matrices

% create all streamline matrices
[ amat, alab ] = feCreateAdjacencyMatrices(pconn, 'all');
[ nmat, nlab ] = feCreateAdjacencyMatrices(pconn, 'nzw');
%[ zmat, zlab ] = feCreateAdjacencyMatrices(pconn, 'nzw_clean');

% combine outputs and labels into 1 matrix
%omat = cat(3, amat, nmat, zmat);
%olab = [ alab nlab zlab ];
omat = cat(3, amat, nmat);
olab = [ alab nlab ];

% network stats
[ glob{1}, node{1}, nets{1} ] = fnNetworkStats(omat(:,:,2), 4);
[ glob{2}, node{2}, nets{2} ] = fnNetworkStats(omat(:,:,10), 4);

%% save whole output

save(output, 'pconn', 'rois', 'omat', 'olab', 'glob', 'node', 'nets', '-v7.3');

% remove parallel pool
delete(pool);

end
