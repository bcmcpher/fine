function [ glob, node, nets ] = fnEstimateLouvainCommunity(mat, gamma, iters, tau)
%fnEstimateLouvainCommunity is a loop to compute and return Louvain
% community estimates with practical iterations in a usable form alongside
% other network summaries.
%
% add small worldness? it's not really a useful measure...
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

% % compute small-worldness
% sw1 = mean(cellfun(@mean, rep.mcoef));
% sw2 = mean(cellfun(@mean, rep.chpl));
% 
% fns.nrm.smwrld = (fns.nrm.mcoef / sw1) / (fns.nrm.chpl / sw2);

end

