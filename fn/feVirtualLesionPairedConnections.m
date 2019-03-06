function [ netw, vltry ] = feVirtualLesionPairedConnections(netw, wlab, M, dsig, nTheta, S0, clobber)
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

if ~isfield(netw.pconn{1}.fibers, wlab)
    error('Fiber weights field ''wlab'' does not exist.');
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

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

% check if vl already exists
if(norm)
    if(isfield(netw.pconn{1}, 'vl_nrm') && clobber == 0)
        error('Virtual lesions on the normalized signal already exist. Please set clobber to 1 to explicitly overwrite existing values.');
    end    
else
    if(isfield(netw.pconn{1}, 'vl_raw') && clobber == 0)
        error('Virtual lesions on the raw signal already exist. Please set clobber to 1 to explicitly overwrite existing values.');
    end
end

% preallocate virtual lesion output
%vlout = cell(length(pconn), 1);
vlcnt = 0;
vltry = 0;

% extract the edge list
pconn = netw.pconn;

disp(['Computing virtual lesion on a possible ' num2str(size(pconn, 1)) ' edges...']);

tic;
for ii = 1:size(pconn, 1)
    
    % pull the edge
    edge = pconn{ii};
    indx = pconn{ii}.fibers.indices;
    wght = pconn{ii}.fibers.(wlab);
    
    % if the weights field is mepty
    if sum(wght) == 0
        
        % set to zero and continue
        edge.vl.s.mean  = 0;
        edge.vl.em.mean = 0;
        edge.vl.j.mean  = 0;
        edge.vl.kl.mean = 0;
                
    else
        
        if(norm)
            
            try
                % compute a normalized virtual lesion
                [ ewVL, ewoVL ] = feComputeVirtualLesionM_norm(M, wght, dsig, nTheta, indx, S0);
                edge.vl         = feComputeEvidence_norm(ewoVL, ewVL);
                vlcnt = vlcnt + 1;
            catch
                vltry = vltry + 1;
                edge.vl.s.mean  = 0;
                edge.vl.em.mean = 0;
                edge.vl.j.mean  = 0;
                edge.vl.kl.mean = 0;
            end
            
        else
            
            try
                % compute a raw virtual lesion
                [ ewVL, ewoVL ] = feComputeVirtualLesionM(M, wght, dsig, nTheta, indx);
                edge.vl         = feComputeEvidence(ewoVL, ewVL);
                vlcnt = vlcnt + 1;
            catch
                vltry = vltry + 1;
                edge.vl.s.mean  = 0;
                edge.vl.em.mean = 0;
                edge.vl.j.mean  = 0;
                edge.vl.kl.mean = 0;
            end
            
        end
        
    end
    
    % reassign edge w/ computed virtual lesion
    pconn{ii} = edge;
    
end
time = toc;

% reassign edges with virtual lesion added
netw.pconn = pconn;

disp(['Computed ' num2str(vlcnt) ' vitual lesions in ' num2str(round(time)/60) ' minutes.']);

clear ii vlcnt ewVL ewoVL time

% % add virtual lesion to pconn object
% for ii = 1:length(pconn)
%     
%     % pull subset field
%     tmp = pconn{ii}.(label);
%     
%     % add vl output as normed or raw
%     if(norm)
%         tmp.vl_nrm = vlout{ii};
%         
%         % assign virtual lesion matrix entries
%         tmp.matrix.soe_nrm = vlout{ii}.s.mean;
%         tmp.matrix.emd_nrm = vlout{ii}.em.mean;
%         tmp.matrix.kl_nrm  = mean(vlout{ii}.kl.mean);
%         tmp.matrix.jd_nrm  = mean(vlout{ii}.j.mean);
%     else
%         tmp.vl_raw = vlout{ii};
%         
%         % assign virtual lesion matrix entries
%         tmp.matrix.soe_raw = vlout{ii}.s.mean;
%         tmp.matrix.emd_raw = vlout{ii}.em.mean;
%         tmp.matrix.kl_raw  = vlout{ii}.kl.mean;
%         tmp.matrix.jd_raw  = vlout{ii}.j.mean;
%     end
%     
%     % reassign virtual lesions to paired connection array
%     pconn{ii}.(label) = tmp;
%     
% end

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
