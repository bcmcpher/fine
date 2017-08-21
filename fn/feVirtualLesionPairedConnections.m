function [ pconn, vlout ] = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, label)
%feVirtualLesionPairedConnections runs virtual lesion on all edges store in pconn.
%
% INPUTS:
%     M                - fe model, sparse tensor object from fit fe structure
%     fascicle_weights - vector of streamline weights from fit fe structure
%     measured_dsig    - measured diffusion signal estimated from fe structure
%     nTheta           - the number of directions stored in an fe structure
%     pconn            - the paired connection object
%     label - string indicating the fiber groups for which to create virtual lesions
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%             Additionally, this can be run after cleaning, resulting in
%             valid calls of 'all_clean' and 'nzw_clean', respectively.
%
% OUTPUTS:
%     pconn - is the paired connections object with the virtual lesion output added
%             to each connection.
%
%     vlout - debugging output; the cell array that is added internally to pconn
%
% TODO:
% - add 'norm' option for computing virtual lesions
% - check that the results I mean() are supposed to be meaned
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
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create a virtual lesion for every connection with non-zero weighted fibers
% pconn = feVirtualLesionPairedConnections(M, fascicle_weights, measured_dsig, nTheta, pconn, 'nzw');
%
% Brent McPherson (c), 2017 - Indiana University
%

% run parallelized virtural lesion

% preallocate virtual lesion output
vlout = cell(length(pconn), 1);
vlcnt = 0;

tic;
parfor ii = 1:length(pconn)
    
    % pull field requested
    tmp = pconn{ii}.(label);
    
    if sum(tmp.weights) == 0
        
        % set to zero and continue
        vlout{ii}.s.mean  = 0;
        vlout{ii}.em.mean = 0;
        vlout{ii}.j.mean  = 0;
        vlout{ii}.kl.mean = 0;
        continue
        
    else
        % CHECK ORDER OF weVL / ewoVL WITH CESAR
        % compute a virtual lesion
        [ ewVL, ewoVL ] = feComputeVirtualLesionM(M, fascicle_weights, measured_dsig, nTheta, tmp.indices);
        vlout{ii}       = feComputeEvidence(ewoVL, ewVL);
        vlcnt = vlcnt + 1;
    end
    
end
time = toc;

display(['Computed ' num2str(vlcnt) ' vitual lesions in ' num2str(round(time)/60) ' minutes.']);

clear ii vlcnt ewVL ewoVL time

% add virtual lesion to pconn object
for ii = 1:length(pconn)
    
    % pull subset field
    tmp = pconn{ii}.(label);
    
    % add whole vl output
    tmp.vl = vlout{ii};
    
    % assign virtual lesion matrix entries
    tmp.matrix.soe = vlout{ii}.s.mean;
    tmp.matrix.emd = vlout{ii}.em.mean;
    tmp.matrix.kl  = mean(vlout{ii}.kl.mean);
    tmp.matrix.jd  = mean(vlout{ii}.j.mean);
    
    % reassign virtual lesions to paired connection array
    pconn{ii}.(label) = tmp;
    
end

end

%% internal fxns to not overflow memory in parpool

function [ rmse_wVL, rmse_woVL, nFib_tract, nFib_PN, nVoxels] = feComputeVirtualLesionM(M, fascicle_weights, measured_dsig, nTheta, ind_tract)
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
ind_nnz   = find(fascicle_weights);
ind2      = intersect(ind_nnz,val);

nFib_tract = length(ind_tract);
nFib_PN    = length(ind2);

%% Compute rmse restricted to Path neighborhood voxels without Virtual Lesion
measured_dsig = measured_dsig(:,voxel_ind);

% Restrict tensor model to the PN voxels
M.Phi   = M.Phi(:,voxel_ind,:);
nVoxels = size(M.Phi,2);

% Generate predicted signal and model error with the full-unleasioned model
w_woVL         = fascicle_weights;
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
w_VL            = fascicle_weights;
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

