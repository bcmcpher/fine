function [ netw ] = fnComputeMatrixField(netw, clobber)
%fnComputeMatrixField creates the valuest that are assigned to adjacency
% matrices for the pconn list object. New measures should be added here.
%
% If pvoxels exists, I can compute volume, ADF, any micro-structure average
% How should this fxn do that? 
% 

%% parse arguments

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

if clobber == 0 && isfield(netw.pconn{1}, 'matrix')
    error('Matrix field has already been computed. Set clobber to true to recompute values.');
end

if clobber == 1 && isfield(netw.pconn{1}, 'matrix')
    disp('Matrix field will be overwritten for every edge...');
end

%% compute the matrix fields

% grab the sorted indice pairs
pairs = netw.parc.pairs;

% for every connection
for ii = 1:size(netw.pconn, 1)
    
    % pull the rois for the pair and the edge information
    roi1 = netw.rois{pairs(ii, 1)};
    roi2 = netw.rois{pairs(ii, 2)};
    edge = netw.pconn{ii};
    
    % calculate combined size of ROI (voxel count / volume in mm^2)
    psz = roi1.size + roi2.size;
    pvl = roi1.volume + roi2.volume;
    
    % pull the streamline count
    cnt = size(edge.fibers.indices, 1);
    
    % pull the average streamline length
    len = mean(edge.fibers.lengths, 'omitnan');
    
    % compute standard error of length
    len_se = std(edge.fibers.lengths, 'omitnan') / sqrt(size(edge.fibers.lengths, 1));
    
    % Replace mean/se of length w/ 0 if it's empty
    if isnan(len)
        len = 0;
        len_se = 0;
    end
        
    % create 1 / sum of lengths for Hagmann's correction
    if isempty(edge.fibers.lengths)
        dln = 0;
    else
        dln = sum(1 / edge.fibers.lengths);
    end

    % create all count edge measures
    matrix.count = cnt;
    matrix.densz = (2 * cnt) / psz;
    matrix.denvl = (2 * cnt) / pvl;
    matrix.length = len;
    matrix.len_se = len_se;
    matrix.densln = (2 / psz) * dln;
    matrix.denvln = (2 / pvl) * dln;
    
    % if norm virtual lesion is present, store the values
    if isfield(edge, 'vl_nrm')
        matrix.soe_nrm_mn = edge.vl_nrm.s.mean;
        matrix.soe_nrm_sd = edge.vl_nrm.s.std;
        matrix.emd_nrm = edge.vl_nrm.em.mean;
        matrix.kld_nrm = edge.vl_nrm.kl.mean;
        matrix.jef_nrm = edge.vl_nrm.j.mean;
    end
    
    % if raw virtual lesion is present, store the values
    if isfield(edge, 'vl_raw')
        matrix.soe_raw_mn = edge.vl_raw.s.mean;
        matrix.soe_raw_sd = edge.vl_raw.s.std;
        matrix.emd_raw = edge.vl_raw.em.mean;
        matrix.kld_raw = edge.vl_raw.kl.mean;
        matrix.jef_raw = edge.vl_raw.j.mean;
    end
    
    % if volume is present, compute the volume measures
    
    % if microstructure is present, compute mn/std microstructural values
    
    % what else is optional?
    
    % assign the output back
    edge.matrix = matrix;
    
    % reassign edge
    netw.pconn{ii} = edge;
    
end

end

