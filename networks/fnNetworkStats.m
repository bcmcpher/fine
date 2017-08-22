function [ glob, node, nets ] = fnNetworkStats(mat)
%fnNetworkStats returns BCT measures of a symetric undirected network.
% Additionally, it computes Small World Propensity from the NCT.
%
% INPUTS:
%     mat     - the input adjacency matrix; can be weighted or unweighted
%
% OUTPUTS:
%     glob - structure conataining global network properties
%     node - structure containing node wise network properties
%     nets - structure containing graphs describing properties
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
% [ glob, node, nets ] = fnNetworkStats(omat(:,:,1));
%
% Brent McPherson (c), 2017 - Indiana University
%

%% fiber data - connection matrix

display('Thresholding Data...');

nets.raw = mat;
nets.nrm = weight_conversion(nets.raw, 'normalize');
nets.len = weight_conversion(nets.raw, 'lengths');
[ nets.dist, nets.edge ] = distance_wei(nets.len); 

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
glob.assort = assortativity_wei(nets.nrm, 0);

% will fail if log10(mat) is provided
[ node.corper, glob.corper ] = core_periphery_dir(nets.nrm, 1); node.corper = node.corper';

[ glob.charpl, glob.efficiency, node.eccentricity, glob.radius, glob.diameter ] = charpath(nets.dist);
glob.rcc = rich_club_wu(nets.nrm);

% compute small world propensity
[ glob.swp, glob.swp_dc, glob.swp_dl, glob.swp_dt ] = rep_swp(nets.raw, 50);

end
