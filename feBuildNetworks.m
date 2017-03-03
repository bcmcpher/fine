function [ emat, cmat, pconn, rois, ematLabs, cmatLabs ] = feBuildNetworks(fe, aparc, msvol, nclust, cacheDir)
%build_networks_demo creates all meaningful combinations of data for my
% current project. I will modify this for future processing as things
% develop.
%
% INPUTS:
% - 'fe' is a fit fe structure that is used to build paired connections indices
% - 'aparc' is the labeled parcellation file used to identify nodes
% - 'msvol' is the aligned microstructural volume to compute tract profiles 
% - 'nclust' defines the number of cores for parpool
% - 'cacheDir' is the explicit tmp directory 
%
% OUTPUTS:
% - 'pconn' is the paired connections object with VL data stored and added
%   to the matrix field for generating networks
% - 'tprof' is the cell array of profiled output only as a cell array. For debugging.
%

%% load data when necessary

% if string is passed, assume it's a path and load it
if isstring(fe) 
    display('Loading fe data...');
    load(fe);
end

% if string is passed, assume it's a path and load it
if isstring(aparc)
    aparc = niftiRead(aparc);
end

% load msvol if it's not already loaded
if isstring(msvol)
    msvol = niftiRead(msvol);
end

%% start parallel pool

% start parallel pool for fefgGet
display(['Opening parallel pool with ', num2str(nclust), ' cores...']);

% create parallel cluster object
c = parcluster;

% set number of cores from arguments
c.NumWorkers = nclust;

% set temporary cache directory
t = tempname(cacheDir);

% make cache dir
OK = mkdir(t);

% check and set cachedir location
if OK
    % set local storage for parpool
    c.JobStorageLocation = t;
end

% start parpool - close parpool at end of fxn
parpool(c, nclust, 'IdleTimeout', 120);

%% build intital paired connection structure

[ pconn, rois ] = feCreatePairedConnections(fe, aparc);

%% perform cleaning on individual connections

% clean all streamline connections
pconn = feCleanPairedConnections(fe, pconn, 'all');

% clean all non-zero weighted streamline connections
pconn = feCleanPairedConnections(fe, pconn, 'nzw');

%% compute virtual lesions on connections

% run virtual lesions on all nzw streamline connections
pconn = feVirtualLesionPairedConnections(fe, pconn, 'nzw');

% run virtual lesions on all cleaned streamline connections
pconn = feVirtualLesionPairedConnections(fe, pconn, 'all_clean');

% run virtual lesions on all cleaned nzw streamline connections
pconn = feVirtualLesionPairedConnections(fe, pconn, 'nzw_clean');

%% compute tract profiles

% estimate FA tract profile for all 
pconn = feTractProfilePairedConnections(fe, pconn, 'all', msvol, 'fa');

% estimate FA tract profile for nzw
pconn = feTractProfilePairedConnections(fe, pconn, 'nzw', msvol, 'fa');

% estimate FA tract profile for all cleaned
pconn = feTractProfilePairedConnections(fe, pconn, 'all_clean', msvol, 'fa');

% estimate FA tract profile for nzw cleaned
pconn = feTractProfilePairedConnections(fe, pconn, 'nzw_clean', msvol, 'fa');

%% create output matrices

[ all, allLabs ] = feCreateAdjacencyMatrices(pconn, 'all');

[ nzw, nzwLabs ] = feCreateAdjacencyMatrices(pconn, 'nzw');

[ alc, alcLabs ] = feCreateAdjacencyMatrices(pconn, 'all_clean');

[ nzc, nzcLabs ] = feCreateAdjacencyMatrices(pconn, 'nzw_clean');

% create output matrix and labels
emat = cat(3, all, nzw);
ematLabs = [ allLabs nzwLabs ];

cmat = cat(3, alc, nzc);
cmatLabs = [ alcLabs nzcLabs ];

end