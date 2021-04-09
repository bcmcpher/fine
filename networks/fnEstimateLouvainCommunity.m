function [ glob, node, nets ] = fnEstimateLouvainCommunity(mat, iters, gamma, tau)
%fnEstimateLouvainCommunity is a loop to compute and return Louvain
% community estimates with a user defined number of iterations alongside
% other network summaries.
%
% INPUTS:
%     mat   - the input adjacency matrix; can be weighted or unweighted
%     iters - the number of iterations to find the maximum Q statistic
%     gamma - Louvain resolution parameter; can be a range of values - the
%             gamma with the maximum q-stat will be returned
%                 gamma > 1:  detects smaller modules
%            0 <= gamma < 1:  detects larger modules
%                 gamma = 1:  Louvain (classic) modularity
%     tau   - threshold of reclustering for consensus
%
% OUTPUTS:
%     glob - structure conataining global network properties of Louvain estimation
%     node - structure containing node wise network properties of Louvain estimation
%     nets - structure containing graphs describing properties of Louvain estimation
%
% TODO:
% - vectorize preallocation of iters / gamma to allow parallel execution
%
% EXAMPLE:
%
% % load data
% parc          = niftiRead('labels.nii.gz');
% fg            = feGet(fe, 'fibers acpc');
% fibers        = fg.fibers;
% fibLength     = fefgGet(fg, 'length');
% weights       = feGet(fe, 'fiberweights');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create adjacency matrices of non-zero weighted fibers
% [ omat, olab ] = feCreateAdjacencyMatrices(pconn, 'nzw');
%
% % compute binary network statistics from omat output
% [ glob, node, nets ] = fnRentianScaling(omat(:,:,1), rois, 50);
%
% Brent McPherson (c), 2017 - Indiana University
%

% parse optional arguments
if(~exist('gamma', 'var') || isempty(gamma))
    gamma = 1;
end

if(~exist('tau', 'var') || isempty(tau))
    tau = 0;
end

% parse gamma for size
ngam = size(gamma, 2);

% create simple data structure
nets.raw = mat;
nets.nrm = weight_conversion(nets.raw, 'normalize');

% perallocate loop
ci = zeros(size(mat, 1), ngam, iters);
q = zeros(ngam, iters, 1);

% for every gamma / iteration
for ii = 1:ngam
    for jj = 1:iters
    
        % estimate community estimates on input data
        [ ci(:, ii, jj), q(ii, jj) ] = community_louvain(mat, gamma(ngam), [], 'modularity');

    end
end

% THIS SHOULD BE THE LONGEST UNCHANGED Q STAT, NOT THE HIGHEST
% take the highest q-value (community structure statistic)
% - find the highest estimated number of neighborhoods
% - sort by neighborhoods

% pull the max q-stat and it's indices in the iterations
% THIS SHOULDN'T BE THE MAX, IT SHOULD BE THE CENTER OF THE MOST FLAT SPOT
[ qmax, vmax ] = max(q(:));
[ g_max, i_max ] = ind2sub(size(q), vmax);

% store the output
glob.qstat = qmax;
glob.gamma = gamma(max(g_max));

% group nodes
node.assign = ci(:, g_max, i_max);
ci_max = squeeze(ci(:, g_max, i_max));

% sort matrix by nodes
[ ~, bb ] = sort(ci_max);

% overlap of nieghborhood assignments
agree = agreement(ci_max);

% normalize for resampling N
agree = agree ./ iters;

% return agreement matrix and sorted normalized data
nets.agree = agree(bb, bb);
nets.rsrta = nets.raw(bb, bb);
nets.nsrta = nets.nrm(bb, bb);

% build consensus partition or agreement matrix
node.consensus = consensus_und(nets.agree, tau, iters);

end

