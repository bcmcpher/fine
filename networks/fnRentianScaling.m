function [ p, se, fh ] = fnRentianScaling(netw, reps)
%fnRentianScaling computes and returns a network measure that controls for
% the physical space of the network in partitioning nodes.
%
% INPUTS:
%     netw - the input network object with values estimated
%     reps - The number of repetitions used in the estimate (default = 5000)
%
% OUTPUTS:
%     p  - Rent's exponent
%     se - standard error estimate from repetitions
%     fh - figure handle of plot of exponent
%
% TODO:
% - figure out the best default range for the plot
% - determine the theoretical minimum given an input?
% - plot the trend line
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

if(~exist('reps', 'var') || isempty(reps))
    reps = 5000;
end

disp('Extracting Network...');

% pull edge weight to be binarized
% this ~should~ always be the same, regardless of edge weight selected
val = cellfun(@(x) length(x.fibers.indices), netw.edges);
mat = squareform(val);

% binarize data matrix
dat = weight_conversion(mat, 'binarize');

disp('Extracting coordinates...');

% grab the roi coordinates
coords = cell2mat(cellfun(@(x) x.center.acpc, netw.nodes, 'UniformOutput', false));

disp('Estimating the Rentian Scaling...');

% rentian scaling - separate fxn
[ N, E ] = rentian_scaling_3d(dat, coords, reps, 0.000001); 
% N is number of nodes in each partition
% E is the number of boundaries crossing each partition

% estimate Rent's exponent
[ b, stats ] = robustfit(log10(N), log10(E));

% return Rent's exponent and standard error from robustfit()
p = b(2);
se = stats.se(2);

% plot from description
fh = figure; 
loglog(E, N, '*');
title([ 'Rent''s Exponent = ' num2str(p) ' +/- ' num2str(se) ' s.e.' ]);

end

