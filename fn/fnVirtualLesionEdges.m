function [ netw ] = fnVirtualLesionEdges(netw, M, weights, dsig, nTheta, S0, nbins, nmc, nboots, thr, clobber)
%feVirtualLesionPairedConnections runs virtual lesion on all edges store in pconn.
%
% INPUTS:
%     pconn         - the paired connection object
%     label         - string indicating the fiber groups for which to create virtual lesions
%                     either:
%                         'all' for all assigned streamlines or
%                         'nzw' for non-zero weighted fibers returned by LiFE
%                     Additionally, this can be run after cleaning, resulting in
%                     valid calls of 'all_clean' and 'nzw_clean', respectively.
%     M             - fe model, sparse tensor object from fit fe structure
%     weights       - vector of streamline weights from fit fe structure
%     measured_dsig - measured diffusion signal estimated from fe structure
%     nTheta        - the number of directions stored in an fe structure
%     S0            - the mean diffusion signal. If this is passed, a
%                     normalized virtual lesion is performed.
%
% OUTPUTS:
%     pconn - is the paired connections object with the virtual lesion output added
%             to each connection.
%
%     vlout - debugging output; the cell array that is added internally to pconn
%
% TODO:
% - add std/pval for SOE to empty conditions
%
% EXAMPLE:
%
% % load data
% parc          = niftiRead('labels.nii.gz');
% fg            = feGet(fe, 'fibers acpc');
% fibers        = fg.fibers;
% fibLength     = fefgGet(fg, 'length');
% weights       = feGet(fe, 'fiberweights');
% nTheta        = feGet(fe, 'nbvals');
% M             = feGet(fe, 'model');
% measured_dsig = feGet(fe, 'dsigdemeaned by voxel');
%
% % if you want the normalized estimate, load the mean dsig
% S0 = feGet(fe, 'b0signalimage');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create a virtual lesion for every connection with non-zero weighted fibers
% pconn = feVirtualLesionPairedConnections(pconn, 'nzw', M, weights, measured_dsig, nTheta);
%
% % create a normalized virtual lesion for every connection with non-zero weighted fibers
% pconn = feVirtualLesionPairedConnections(pconn, 'nzw', M, weights, measured_dsig, nTheta, S0);
%
% Brent McPherson & Franco Pestilli (c), 2017 - Indiana University
%

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = false;
end

% check if vl already exists, requires clobber to be set to true to rerun
if(isfield(netw.parc, 'vl') && ~clobber)
    error('Virtual lesions already exist. Please set clobber to 1 to explicitly overwrite existing values.');
end

% parse optional arguments to determine if regular or normed vl is run
if(~exist('S0', 'var') || isempty(S0))
    S0 = [];
    norm = 0;
    netw.vl.sig = 'raw';
    disp('Computing virtual lesions on raw diffusion signal.');
else
    norm = 1;
    netw.vl.sig = 'normalized';
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

% store vl parameters in network object
netw.vl.nbin = nbins;
netw.vl.nmc = nmc;
netw.vl.nboots = nboots;
netw.vl.thr = thr;

% preallocate virtual lesion output
vlcnt = 0;
vltry = 0;

% extract the edge list
pconn = netw.edges;

% I need a better estimate of time - use this for now
if isfield(pconn{1}.fibers, 'weights')
    wconn = cellfun(@(x) any(x.fibers.weights > 0), pconn);
    disp(['Computing ' num2str(sum(wconn)) ' virtual lesions of a possible ' num2str(size(pconn, 1)) ' edges...']);
else
    disp(['Computing virtual lesions on a possible ' num2str(size(pconn, 1)) ' edges...']);
end

tic;
for ii = 1:size(pconn, 1)
    
    % pull the edge and weights from passed vector
    % - could rely on stored weights, but likely to cause confusion
    edge = pconn{ii};
    wght = weights(edge.fibers.indices);
    
    % if there is not any weight greater than 0
    if ~any(wght > 0)
        
        % set to zero and continue
        edge.vl.s.mean  = 0;
        edge.vl.s.std   = 0;
        edge.vl.s.tmean = 0;
        edge.vl.em.mean = 0;
        edge.vl.j.mean  = 0;
        edge.vl.kl.mean = 0;
        
    else
        
        % only keep weighted indices
        keep = wght > 0;
        indx = edge.fibers.indices(keep);
        
        % safely run b/c this is touchy at times
        try
            if(norm)                
                % compute a normalized virtual lesion
                [ ewVL, ewoVL ] = feComputeVirtualLesionM_norm(M, weights, dsig, nTheta, indx, S0);
            else
                % compute the raw virtual lesion
                [ ewVL, ewoVL ] = feComputeVirtualLesionM(M, weights, dsig, nTheta, indx);
            end
            
            % compute evidence based on virtual lesion
            edge.vl = feComputeEvidence_new(ewoVL, ewVL, nbins, nmc, nboots, thr);
            vlcnt = vlcnt + 1;
            
        catch
            
            % add to try count and fill in nan not 0
            vltry = vltry + 1;
            edge.vl.s.mean  = nan;
            edge.vl.s.std   = nan;
            edge.vl.s.tmean = nan;
            edge.vl.em.mean = nan;
            edge.vl.j.mean  = nan;
            edge.vl.kl.mean = nan;
            
        end
        
    end
        
    % reassign edge w/ computed virtual lesion
    pconn{ii} = edge;
    
end
time = toc;

% reassign edges with virtual lesion added
netw.edges = pconn;

disp(['Computed ' num2str(vlcnt) ' vitual lesions in ' num2str(round(time)/60) ' minutes.']);

if vltry > 0
    disp(['Attempted and failed to estimate ' num2str(vltry) ' vitual lesions.']);
end

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
