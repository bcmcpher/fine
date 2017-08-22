function [ p, se, fh ] = fnRentianScaling(mat, rois, reps)
%fnRentianScaling computes and returns a network measure that controls for
% the physical space of the network in partitioning nodes.
%
% INPUTS:
%     mat  - the input adjacency matrix; can be weighted or unweighted
%     rois - ROI structure containing coordinate centers of nodes for scaling
%     reps - The number of repetitions used in the estimate
%
% OUTPUTS:
%     p  - Rent's exponent
%     se - standard error estimate from repetitions
%     fh - figure handle of plot of exponent
%
% TODO:
% - set up a defualt number of iterations
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

% binarize input
dat = weight_conversion(mat, 'binarize');

% grab node coords from rois object 
coords = zeros(length(rois), 3);
for ii = 1:length(rois)
    coords(ii, :) = rois{ii}.centroid.acpc;
end

% rentian scaling - separate fxn
[ N, E ] = rentian_scaling_3d(dat, coords, reps, 0.000001); % n is number of partitions

% estimate Rent's exponent
[ b, stats ] = robustfit(log10(N), log10(E));

% return Rent's exponent and standard error from robustfit()
p = b(2);
se = stats.se(2);

% plot from description
fh = figure; 
loglog(E, N, '*');

end

