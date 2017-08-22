function [ glob, node, nets ] = fnEstimateLouvainCommunity(mat, iters, gamma, tau)
%fnEstimateLouvainCommunity is a loop to compute and return Louvain
% community estimates with a user defined number of iterations alongside
% other network summaries.
%
% INPUTS:
%     mat   - the input adjacency matrix; can be weighted or unweighted
%     iters - the number of iterations to find the maximum Q statistic
%     gamma - Louvain resolution parameter
%                 gamma > 1:  detects smaller modules
%            0 <= gamma < 1:  detects larger modules
%                 gamma = 1:  classic modularity
%     tau   - threshold of reclustering for consensus
%
% OUTPUTS:
%     glob - structure conataining global network properties of Louvain estimation
%     node - structure containing node wise network properties of Louvain estimation
%     nets - structure containing graphs describing properties of Louvain estimation
%
% TODO:
% - add small worldness estimate in with repeats
% - - make this optional? very slow in large networks, extra measure might not be useful
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

% create simple data structure
nets.raw = mat;
nets.nrm = weight_conversion(nets.raw, 'normalize');

% perallocate loop
ci = zeros(size(mat, 1), iters);
q = zeros(iters, 1);

% for every iteration
for ii = 1:iters
    
    % estimate community estimates on input data
    [ ci(:, ii), q(ii) ] = community_louvain(mat, gamma, [], 'modularity');

end;

% take the highest q-value (community structure statistic)
% - find the highest estimated number of neighborhoods
% - sort by neighborhoods
[ ~, indx_max ] = max(q);
glob.qstat = max(q);

% group nodes
ci_max = ci(:, indx_max);
node.assign = ci_max;

% sort matrix by nodes
[ ~, bb ] = sort(ci_max);

% overlap of nieghborhood assignments
agree = agreement(ci);

% normalize for resampling N
agree = agree ./ iters;

% return agreement matrix and sorted normalized data
nets.agree = agree(bb, bb);
nets.rsrta = nets.raw(bb, bb);
nets.nsrta = nets.nrm(bb, bb);

% build consensus partition or agreement matrix
node.consensus = consensus_und(nets.agree, tau, iters);

end

