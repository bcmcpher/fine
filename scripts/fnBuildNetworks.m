function [ pconn, rois, omat, olab ] = fnBuildNetworks(fe, labels, fa, output, nclust, cacheDir)
% fnBuildNetworks creates all the network data from any fit fe and parcellation.
%
% INPUTS: all inputs are strings
%
%   fe       - string containing path name to .mat containing only a fit fe structure
%   labels   - string containing path to nifti of labeled ROIs
%   fa       - string containing path to fa volume for tract profiles
%   output   - string containing path to .mat to save network data
%   nclust   - number of cores for parpool
%   cacheDir - where parpool can be opened from qsub
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

% load labeled aparc+aseg volume and fa map
parc = niftiRead(labels);

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

% Run operations
[ pconn, rois ] = feCreatePairedConnections(parc, fg.fibers, fascicle_length, fascicle_weights);

% virtual lesion matrix
pconn = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, 'nzw');

% compute tract profiles
pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', favol, 'fa');

%% create adjacency matrices

% create all streamline matrices
[ amat, alab ] = feCreateAdjacencyMatrices(pconn, 'all');
[ nmat, nlab ] = feCreateAdjacencyMatrices(pconn, 'nzw');

% combine outputs and labels into 1 matrix
omat = cat(3, amat, nmat);
olab = [ alab nlab ];

% add stats 
[ glob{1}, node{1}, nets{1} ] = fnNetworkStats(omat(:,:,2));
[ glob{2}, node{2}, nets{2} ] = fnNetworkStats(omat(:,:,10));

%% save whole output

% write output to .mat
save(output, 'pconn', 'omat', 'olab', 'glob', 'node', 'nets', '-v7.3');

% remove parallel pool
delete(pool);

end
