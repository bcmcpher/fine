function [ glob, node, nets ] = fnNetworkStats(mat)
%fnNetworkStats returns BCT measures of a symetric undirected network
%
% thresholding / log transform should be preformed outside this fxn  
% should normalization be expected before as well?
%

%% fiber data - connection matrix

display('Thresholding Data...');

nets.raw = mat;
nets.nrm = weight_conversion(nets.raw, 'normalize');
nets.len = weight_conversion(nets.raw, 'lengths');
[ nets.dist, nets.edge ] = distance_wei(nets.len); 
%[ nets.tree, nets.clus ] = backbone_wu(nets.nrm, score);

%% node measures

% compute node network statistics
display('Computing Nodal Summary Statistics...');

node.degree = degrees_und(nets.nrm)';
node.strength = strengths_und(nets.nrm)';
node.between = betweenness_wei(nets.dist);
[ nets.btwcm, node.btwcv ] = edge_betweenness_wei(nets.dist);
node.locEff = efficiency_wei(nets.nrm, 1);
node.ccoef = clustering_coef_wu(nets.nrm);
node.eigv = eigenvector_centrality_und(nets.nrm); % will fail if log10(mat) is provided
node.pagerank = pagerank_centrality(nets.nrm, 0.85, ones(size(nets.nrm, 1), 1)); % is this reasonable?

%% global measures

% compute global network statistics
display('Computing Global Summary Statistics...');

[ glob.density, glob.dnvrt, glob.dnedg ] = density_und(nets.nrm);
glob.mean_ccoef = mean(node.ccoef);
glob.trans = transitivity_wu(nets.nrm);
glob.glbEff = efficiency_wei(nets.nrm, 0);
%[ nets.score, glob.score ] = score_wu(nets.nrm, score); 
glob.assort = assortativity_wei(nets.nrm, 0);

% will fail if log10(mat) is provided
[ node.corper, glob.corper ] = core_periphery_dir(nets.nrm, 1); node.corper = node.corper';

[ glob.charpl, glob.efficiency, node.eccentricity, glob.radius, glob.diameter ] = charpath(nets.dist);
glob.rcc = rich_club_wu(nets.nrm);

% compute small world propensity
[ glob.swp, glob.swp_dc, glob.swp_dl, glob.swp_dt ] = rep_swp(nets.raw, 50);

%% estimate community structure - compartmentalized
% This iterates, so it's slow in large networks. Proabably keep it separate.
% 
% display('Estimating Louvain Community Structure...');
% 
% [ lglob, lnode, lnets ] = fnEstimateLouvainCommunity(nets.nrm, iters, 0.5);
% 
% % catch data in return structures
% glob.qstat = lglob.qstat;
% node.assign = lnode.assign;
% node.consensus = lnode.consensus;
% nets.agree = lnets.agree;
% nets.rsrta = lnets.rsrta;
% nets.nsrta = lnets.nsrta;
% 
% % community structure - local / global measures
% [ node.modularity, glob.modulartiy ] = modularity_und(nets.nrm, 1);
% node.pcoef = participation_coef(nets.nrm, node.assign);
% node.modz = module_degree_zscore(nets.nrm, node.assign);

end
