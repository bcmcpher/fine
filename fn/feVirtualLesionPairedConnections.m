function [ pconn, vlout ] = feVirtualLesionPairedConnections(pconn, label, M, weights, measured_dsig, nTheta, S0, clobber)
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
% - none
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
% Brent McPherson (c), 2017 - Indiana University
%

% parse optional arguments to determine if regular or normed vl is run
if(~exist('S0', 'var') || isempty(S0))
    S0 = [];
    norm = 0;
    disp('Computing virtual lesions on raw diffusion signal...');
else
    norm = 1;
    disp('Computing virtual lesions on demeaned diffusion signal...');
end

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

% check if vl already exists
if(norm)
    if (isfield(pconn{1}.(label), 'vl_nrm') && clobber == 0)
        error('Normed virtual lesions of these parameters already exist. Please set clobber to 1 to explicitly overwrite existing values.');
    end
    
else
    if (isfield(pconn{1}.(label), 'vl_raw') && clobber == 0)
        error('Raw virtual lesions of these parameters already exist. Please set clobber to 1 to explicitly overwrite existing values.');
    end
end

% preallocate virtual lesion output
vlout = cell(length(pconn), 1);
vlcnt = 0;

tic;
parfor ii = 1:length(pconn)
    
    % pull field requested
    tmp = pconn{ii}.(label);
    
    if sum(tmp.weights) == 0
        
        % set to zero and continue
        vlout{ii}.s.mean = nan;
        vlout{ii}.s.std = nan;
        vlout{ii}.s.pval = nan;
        vlout{ii}.em.mean = nan;
        vlout{ii}.em.se = nan;
        vlout{ii}.em.pval = nan;
        vlout{ii}.j.mean = nan;
        vlout{ii}.j.se = nan;
        vlout{ii}.j.pval = nan;
        vlout{ii}.kl.mean = nan;
        vlout{ii}.kl.se = nan;
        vlout{ii}.kl.pval = nan;
        vlout{ii}.dis.mean = nan;
        vlout{ii}.dis.se = nan;
        vlout{ii}.dis.pval = nan;
        vlout{ii}.cos.mean = nan;
        vlout{ii}.cos.se = nan;
        vlout{ii}.cos.pval = nan;
        
        continue
        
    else
        
        if(norm)
            
            % compute a normalized virtual lesion
            [ ewVL, ewoVL ] = feComputeVirtualLesionM_norm(M, weights, measured_dsig, nTheta, tmp.indices, S0);
            vlout{ii}       = feComputeEvidence_new(ewoVL, ewVL);
            vlcnt = vlcnt + 1;
            
        else
 
            % compute a raw virtual lesion
            [ ewVL, ewoVL ] = feComputeVirtualLesionM(M, weights, measured_dsig, nTheta, tmp.indices);
            vlout{ii}       = feComputeEvidence(ewoVL, ewVL); % needs to be updated to match _norm fxn
            vlcnt = vlcnt + 1;
        end
        
    end
    
end
time = toc;

disp(['Computed ' num2str(vlcnt) ' vitual lesions in ' num2str(round(time)/60) ' minutes.']);

clear ii vlcnt ewVL ewoVL time

% add virtual lesion to pconn object
for ii = 1:length(pconn)
    
    % pull subset field
    tmp = pconn{ii}.(label);
    
    % add vl output as normed or raw
    if(norm)
        tmp.vl_nrm = vlout{ii};
        
        % assign virtual lesion matrix entries
        tmp.matrix.soe_nrm = vlout{ii}.s.mean;
        tmp.matrix.std_nrm = vlout{ii}.s.std;
        tmp.matrix.spv_nrm = vlout{ii}.s.pval;
        
        tmp.matrix.emd_nrm = vlout{ii}.em.mean;
        tmp.matrix.kl_nrm  = mean(vlout{ii}.kl.mean);
        tmp.matrix.jd_nrm  = mean(vlout{ii}.j.mean);
        
        tmp.matrix.dis_nrm = vlout{ii}.dis.mean;
        tmp.matrix.cos_nrm = vlout{ii}.cos.mean;
        
    else
        tmp.vl_raw = vlout{ii};
        
        % assign virtual lesion matrix entries
        tmp.matrix.soe_raw = vlout{ii}.s.mean;
        
        tmp.matrix.emd_raw = vlout{ii}.em.mean;
        tmp.matrix.kl_raw  = vlout{ii}.kl.mean;
        tmp.matrix.jd_raw  = vlout{ii}.j.mean;
        
        % not present in not feComputeEvidence yet
        %tmp.matrix.dis_raw = vlout{ii}.dis.mean;
        %tmp.matrix.cos_raw = vlout{ii}.cos.mean;
        
    end
    
    % reassign virtual lesions to paired connection array
    pconn{ii}.(label) = tmp;
    
end

end

%% internal fxns to not overflow memory in parpool

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

%% Compute rmse restricted to Path neighborhood voxels without Virtual Lesion
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

%% Compute rmse restricted to Path neighborhood voxels without Virtual Lesion
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
