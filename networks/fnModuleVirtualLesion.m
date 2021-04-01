function [ vmat, olab, indx, pcmodl ] = fnModuleVirtualLesion(netw, srt, labs, M, weights, dsig, nTheta, S0, nbins, nmc, nboots, thr)
%fnModuleVirtualLesion() perfurms virtual lesion on all connections within
% or between a network module(s).
%   
% INPUTS:
%     fe    - a fit LiFE model
%     pconn - a paired connection structure assigning every streamline 
%             to a connection
%     label - which set of indices to assign from pconn. Either
%               - all: all streamlines assigned
%               - nzw: non-zero weighted
%               - zwr: zero weighted / removed (cannot virtual lesion)
%     srt   - node community assignment, size(pconn, 1) x 1 array
%             determines the module community of every connection
%     norm  - (optional, default = 0) 1 or 0, whether or not to normalize 
%             the virtual lesion computation
%
% OUTPUTS:
%    vmat   - a nmod x nmod x 4 array of the requested modules with each
%             virtual lesion measure 
%    olab   - a 1 x 4 string array indicating the measure along the 4th
%             dimension of vmat
%    indx   - the raw data used within the function to compute the vl
%    pcmodl - an array with the module assignment of every connection in pconn
%
% Brent McPherson, (c) 2017, Indiana University
%
% TODO
% - can the logic be cleaned up to be more clear
% - more useful printing updates out
%

% are these interesting / necessary?
if(~exist('labs', 'var') || isempty(labs))
    labs = cell(size(netw.node));    
end

% parse optional arguments to determine if regular or normed vl is run
if(~exist('S0', 'var') || isempty(S0))
    S0 = [];
    norm = 0;
    disp('Computing virtual lesions on raw diffusion signal.');
else
    norm = 1;
    disp('Computing virtual lesions on demeaned diffusion signal.');
end

if(~exist('nbin', 'var') || isempty(nbins))
    nbins = 128;    
end

if(~exist('nmc', 'var') || isempty(nmc))
    nmc = 5;
end

if(~exist('nboots', 'var') || isempty(nboots))
    nboots = 250; 
end

if(~exist('thr', 'var') || isempty(thr))
    thr = 0.05; 
end

%% begin pulling data

% compute the size of node assignments and the number of connections
szsrt = length(srt);
nnodes = size(netw.nodes, 1);

% error if the wrong number of labels is passed
if szsrt ~= nnodes
    error('The number of nodes assigned to modules does not match the number of nodes in the network.');
end

% pull the edge pairs object
pcindx = netw.parc.pairs;

% size of pairs object
szpci = size(pcindx, 1);

% find maximum number of modules
nmod = size(unique(srt), 1);

% find upper diagonal (between modules)
uprd = nchoosek(1:nmod, 2);

% find diagonal (within modules)
diag = [ 1:nmod; 1:nmod ]';

% combine the module indices to compare
pairs = [ diag; uprd ];

clear diag uprd

% size of unique module combinations
szprs = size(pairs, 1);

disp(['Assigning all connections to one of ' num2str(szprs) ' unique modules...']);

% preallocate connection assignment array 
pcmodl = nan(szpci, 1);

% catch labels corresponding to module index
indx = repmat({struct('mod1', [], 'mod2', [], 'indices', [])}, [szprs, 1]);

% grab module assignment of every connection
for conn = 1:szpci
    
    % grab the module labels - srt value
    roi1 = srt(pcindx(conn, 1));
    roi2 = srt(pcindx(conn, 2));
    
    % grab and return the module this connection is in
    % this will try both orders of ROI indexing (only 1 will work)
    try
        % use the row index in pairs to assign the community label
        [ ~, pcmodl(conn) ] = intersect(pairs, [ roi1 roi2 ], 'rows');
    catch
        [ ~, pcmodl(conn) ] = intersect(pairs, [ roi2 roi1 ], 'rows');
    end
    
    % grab module index
    mod = pcmodl(conn);
    
    % assign module pairs - is reassigned every loop...
    indx{mod}.mod1 = pairs(mod, 1);
    indx{mod}.mod2 = pairs(mod, 2);
    
    % grab the indices / lengths / weights for every streamline in each connection
    indx{mod}.indices = [ indx{mod}.indices; netw.edges{conn}.fibers.indices ];
    
end

clear conn roi1 roi2 mod

% total summaries
tcon = 0;
tfib = zeros(szprs, 1);

% report the number of connections and total streamlines of each module
for mod = 1:szprs
    
    % pull values from indx
    nconn = size(find(pcmodl == mod), 1);
    nfib = size(indx{mod}.indices, 1);
    
    % print out each module
    disp([ 'Module ' num2str(mod) ' has ' num2str(nconn) ' connections and ' num2str(nfib) ' streamlines.' ]);
    
    % iterate total counts
    tcon = tcon + nconn;
    tfib(mod) = nfib; 
    
end

% report total connections and streamlines assigned (should be all of them)
disp([ 'A total of ' num2str(tcon) ' connections and ' num2str(sum(tfib)) ' streamlines were assigned.' ]);  

clear tcon tfib mod nconn nfib

disp(['Computing virtual lesion for every one of ' num2str(size(indx, 1)) ' network modules...']);

for vl = 1:szprs
    
    if all(weights(indx{vl}.indices) == 0)
        
        warning('Connection %d has no positive weighted streamlines. VL not computed.', vl)
        indx{vl}.vl.em.mean = 0;
        indx{vl}.vl.s.mean = 0;
        indx{vl}.vl.kl.mean = 0;
        indx{vl}.vl.j.mean = 0;
        
        continue
    
    end
        
    % compute virtual lesion as either raw or normalized
    switch norm
        case 0
            [ wVL, woVL ] = feComputeVirtualLesionM(M, weights, dsig, nTheta, indx{vl}.indices);
            indx{vl}.vl = feComputeEvidence_new(woVL, wVL, nbins, nmc, nboots, thr);
        case 1
            [ wVL, woVL ] = feComputeVirtualLesionM_norm(M, weights, dsig, nTheta, indx{vl}.indices, S0);
            indx{vl}.vl = feComputeEvidence_new(woVL, wVL, nbins, nmc, nboots, thr);
        otherwise
            error('Bad VL normalization call that can''t happen?');
    end
    
end

% preallocate output
vmat = zeros(nmod, nmod, 4);

% collect virtual lesions into output matrix
for vl = 1:szprs
    
    % create upper diagonal output matrix
    vmat(pairs(vl, 1), pairs(vl, 2), 1) = indx{vl}.vl.em.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 2) = indx{vl}.vl.s.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 3) = indx{vl}.vl.kl.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 4) = indx{vl}.vl.j.mean;

    % create lower diagonal output matrix
    vmat(pairs(vl, 2), pairs(vl, 1), 1) = indx{vl}.vl.em.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 2) = indx{vl}.vl.s.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 3) = indx{vl}.vl.kl.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 4) = indx{vl}.vl.j.mean;

end

% output labels that do not change
olab = {'EMD', 'SOE', 'KLD', 'JfD'};

end

%% internal fxns to not overflow memory

function [ rmse_wVL, rmse_woVL, nFib_tract, nFib_PN, nVoxels] = feComputeVirtualLesionM(M, weights, measured_dsig, nTheta, ind_tract)
% This function compute the rmse in a path neighborhood voxels with and
% without Virtual Lesion
%
% This function with name ending in M attempts to use only the sparse tensor as input not the full FE strcuture. 
%
% INPUTS:
% fe: fe structure
% ind1: indices to fibers in the tract to be virtually lesioned
%
%  Copyright (2016), Franco Pestilli (Indiana University) - Cesar F. Caiafa (CONICET)
%  email: frakkopesto@gmail.com and ccaiafa@gmail.com

% ind_nnz = find(fe.life.fit.weights);
% ind_tract = ind_nnz(ind1);

% We want to find which voxels that a group of fibers crosses.
[inds, ~] = find(M.Phi(:,:,ind_tract)); % find nnz entries of subtensor
voxel_ind = unique(inds(:,2)); clear inds
[inds, ~] = find(M.Phi(:,voxel_ind,:));
val       = unique(inds(:,3));
val       = setdiff(val,ind_tract);
ind_nnz   = find(weights);
ind2      = intersect(ind_nnz,val);

nFib_tract = length(ind_tract);
nFib_PN    = length(ind2);

% Compute rmse restricted to Path neighborhood voxels without Virtual Lesion
measured_dsig = measured_dsig(:,voxel_ind);

% Restrict tensor model to the PN voxels
M.Phi   = M.Phi(:,voxel_ind,:);
nVoxels = size(M.Phi,2);

% Generate predicted signal and model error with the full-unleasioned model
w_woVL         = weights;
Mw             =  M_times_w(M.Phi.subs(:,1), ...
                            M.Phi.subs(:,2), ...
                            M.Phi.subs(:,3), ...
                            M.Phi.vals, ...
                            M.DictSig, ...
                            w_woVL, ...
                            nTheta, ...
                            nVoxels);
predicted_woVL = reshape(Mw, size(measured_dsig));
rmse_woVL      = sqrt(mean(( measured_dsig - predicted_woVL ).^2 ,1));

% Generate predicted signal and model error with the leasioned model
w_VL            = weights;
w_VL(ind_tract) = 0;
Mw = M_times_w(M.Phi.subs(:,1), ...
               M.Phi.subs(:,2), ...
               M.Phi.subs(:,3), ...
               M.Phi.vals, ...
               M.DictSig, ...
               w_VL, ...
               nTheta, ...
               nVoxels);
predicted_VL = reshape(Mw,size(measured_dsig));
rmse_wVL     = sqrt(mean(( measured_dsig - predicted_VL ).^2 ,1));

end

function [ rmse_wVL, rmse_woVL, nFib_tract, nFib_PN, nVoxels] = feComputeVirtualLesionM_norm(M, weights, measured_dsig, nTheta, ind_tract, S0)
% This function compute the rmse in a path neighborhood voxels with and
% without Virtual Lesion while correcting for the mean diffusion signal
%
% This function with name ending in M attempts to use only the sparse tensor as input not the full FE strcuture. 
%
% INPUTS:
% fe: fe structure
% ind1: indices to fibers in the tract to be virtually lesioned
%
%  Copyright (2016), Franco Pestilli (Indiana University) - Cesar F. Caiafa (CONICET)
%  email: frakkopesto@gmail.com and ccaiafa@gmail.com

% ind_nnz = find(fe.life.fit.weights);
% ind_tract = ind_nnz(ind1);

% We want to find which voxels that a group of fibers crosses.
[inds, ~]  = find(M.Phi(:,:,ind_tract)); % find nnz entries of subtensor
voxel_ind  = unique(inds(:,2)); clear inds
[inds, ~]  = find(M.Phi(:,voxel_ind,:));
val        = unique(inds(:,3));
val        = setdiff(val,ind_tract);
ind_nnz    = find(weights);
ind2       = intersect(ind_nnz,val);

S0 = S0(voxel_ind);

nFib_tract = length(ind_tract);
nFib_PN    = length(ind2);

% Compute rmse restricted to Path neighborhood voxels without Virtual Lesion
measured_dsig = measured_dsig(:,voxel_ind);

% Restrict tensor model to the PN voxels
M.Phi   = M.Phi(:,voxel_ind,:);
nVoxels = size(M.Phi,2);

% Generate predicted signal and model error with the full-unleasioned model
w_woVL         = weights;
Mw             =  M_times_w(M.Phi.subs(:,1), ...
                            M.Phi.subs(:,2), ...
                            M.Phi.subs(:,3), ...
                            M.Phi.vals, ...
                            M.DictSig, ...
                            w_woVL, ...
                            nTheta, ...
                            nVoxels);
predicted_woVL = reshape(Mw, size(measured_dsig));
rmse_woVL      = sqrt(mean(( measured_dsig - predicted_woVL ).^2 ,1));
rmse_woVL      = rmse_woVL ./ S0';

% Generate predicted signal and model error with the leasioned model
w_VL            = weights;
w_VL(ind_tract) = 0;
Mw = M_times_w(M.Phi.subs(:,1), ...
               M.Phi.subs(:,2), ...
               M.Phi.subs(:,3), ...
               M.Phi.vals, ...
               M.DictSig, ...
               w_VL, ...
               nTheta, ...
               nVoxels);
predicted_VL = reshape(Mw,size(measured_dsig));
rmse_wVL     = sqrt(mean(( measured_dsig - predicted_VL ).^2 ,1));
rmse_wVL     = rmse_wVL ./ S0';

end
