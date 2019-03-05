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
    
    % pull the average streamline length, replace w/ 0 if it's empty
    len = mean(edge.fibers.lengths, 'omitnan');
    if isnan(len)
        len = 0;
    end
    
    % create 1 / sum of lengths for Hagmann's correction
    if isempty(edge.fibers.lengths)
        dln = 0;
    else
        dln = sum(1 / edge.fibers.lengths);
    end

    % create all edge measures
    matrix.count = cnt;
    matrix.densz = (2 * cnt) / psz;
    matrix.denvl = (2 * cnt) / pvl;
    matrix.length = len;
    matrix.densln = (2 / psz) * dln;
    matrix.denvln = (2 / pvl) * dln;
            
    % if volume is present, compute the volume measure
    
    % if microstructure is present, compute mn/std microstructural values
    
    % what else is optional?
    
    % assign the output back
    edge.matrix = matrix;
    
    % reassign edge
    netw.pconn{ii} = edge;
    
end

end

