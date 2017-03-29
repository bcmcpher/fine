function [ glob, node, nets ] = fnBinaryNetworkStats(mat, kcore, cliques)
%fnBinaryNetworkStats returns BCT measures of a symetric binary network.
% Because I am always working with weighted networks, this fxn should only
% focus on binary measures that may be useful / desired without repeating
% network measures that should probably be estimated in the directed
% network. 
%
% I will thoroughly note any exceptions.
%
% thresholding / log transform should be preformed outside this fxn  

%% fiber data - connection matrix

display('Thresholding Data...');

nets.raw = mat;
nets.bin = weight_conversion(nets.raw, 'binarize');
nets.len = weight_conversion(nets.bin, 'lengths');
[ nets.dist, nets.edge ] = distance_wei(nets.len); 
[ nets.tree, nets.clus ] = backbone_wu(nets.bin, kcore);

%% compute network statistics

display('Computing Binary Network Summary Statistics...');
 
% stats on binarized matrix / k-core
[ nets.kcore, glob.kn, glob.plord, glob.pllvl ] = kcore_bu(nets.bin, kcore);
[ node.core, node.ckn ] = kcoreness_centrality_bu(nets.bin);

% reorient kcore stats
node.core = node.core';
node.ckn  = node.ckn';

% generalized topological measure
nets.gtom = gtom(nets.bin, 1);
nets.gdst = 1 - nets.gtom;

% neighbor overlap
[ nets.novr, node.nec, node.ndeg ] = edge_nei_overlap_bu(nets.bin);

% matching index
nets.match = matching_ind_und(nets.bin);

% clique communities
node.cliques = clique_communities(nets.bin, cliques); node.cliques = node.cliques';

% subgraph_centrality
node.subgraph = subgraph_centrality(nets.bin);

% joint degree
[ nets.joint, glob.jdOdd, glob.jdIn, glob.jdOut ] = jdegree(nets.bin);

end
