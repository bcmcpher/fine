function [ omat, out ] = fnCreateLinkNetwork(pconn, label)
%fnCreateLinkNetwork creates a link network from a pconn list that has had
% edge volumes precomputed for the label requested.
%
% INPUTS:
%     pconn - is the paired connections object to create a link network from. 
%             Edge volumes need to be estimated beforehand with fnFindPathVoxels.m 
%
%     label - string indicating the fiber groups for which to create virtual lesions
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%             Additionally, this can be run after cleaning, resulting in
%             valid calls of 'all_clean' and 'nzw_clean', respectively.
%
% OUTPUTS:
%     omat - a matrix that is (length(pconn) x length(pconn)) in size
%            showing the Dice coefficient of 
%     out  - cell array of the data used to compute edge entry in omat
%
% TODO:
% - add other metrics?
% Mutual Information and Joint Entropy
% Cosine Similarity - pdist(ld, 'cosine'); or pdist2(ld1, ld2, 'cosine'); 
% where ld is matrix of node indices, counts, and/or dictionary orientations.
% some of these? all of these?
% https://en.wikipedia.org/wiki/Diversity_index#Simpson_index
% Simpson's Index: probability that 2 random voxels are part of intersection
% Richness: how many links does an intersection contain
% Shannon Diversity Index: global measure of intersections between links 
%
% EXAMPLE:
%
% % load data
% parc          = niftiRead('labels.nii.gz');
% fg            = feGet(fe, 'fibers acpc');
% fibers        = fg.fibers;
% fibLength     = fefgGet(fg, 'length');
% weights       = feGet(fe, 'fiberweights');
% Phi           = feGet(fe, 'Phi');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % find the voxels that each edge occupies
% pconn = fnFindPathVoxels(Phi, pconn, 'nzw');
%
% % create link network based on overlap of non-zero weighted fibers
% [omat, out ] = fnCreateLinkNetwork(pconn, 'nzw');
%
% Brent McPherson (c), 2017 - Indiana University
%

%% generate link network

display('Preallocating output structures...');

% number of edges; link nodes
lnodes = size(pconn, 1);

% create indices for every pair of link connections
pairs = nchoosek(1:lnodes, 2);

% build empty cells to parallelize computation
out = cell(length(pairs), 1);

% build empty matrices
omat = zeros(lnodes, lnodes);

display('Computing edge intersections...');

% for every connections overlap
tic;
parfor ii = 1:length(pairs)
    
    % simple index for each edge combination
    ti1 = pairs(ii, 1);
    ti2 = pairs(ii, 2);
    
    % pull labels
    out{ii}.edge1 = ti1;
    out{ii}.edge2 = ti2;
    
    % grab the unique voxels of each edge
    li1 = pconn{ti1}.(label).pvoxels;
    li2 = pconn{ti2}.(label).pvoxels;
    
    % if either connection is empty, fill in 0 and move on
    if isempty(li1) || isempty(li2)
        out{ii}.indices = [];
        out{ii}.dice = 0;
        continue
    end
    
    % find the link intersection
    indices = intersect(li1, li2);
    
    % if the intersection is empty, move on
    if isempty(indices)
        out{ii}.indices = [];
        out{ii}.dice = 0;
        continue
    end
    
    % find the size of the intersection
    num = size(indices, 1);
    
    % combine the voxel indices for denominator of Dice coeff
    den = size(li1, 1) + size(li2, 1);
    
    % build Dice coeff values for assignment
    val = (2 * num) / den;
    
    % catch output
    out{ii}.indices = indices;
    out{ii}.dice = val;
    
end
time = toc;

display(['Created link network from ' num2str(length(pairs)) ' unique edge combinations in ' num2str(round(time)/60) ' minutes.']);

%% assemble matrix

display('Creating link network matrix...');

for ii = 1:length(pairs)
    
    % simple index for each 
    ti1 = pairs(ii, 1);
    ti2 = pairs(ii, 2);
    
    % build matrix from out
    omat(ti1, ti2) = out{ii}.dice;
    omat(ti2, ti1) = out{ii}.dice;
    
end

clear ii ti1 ti2 

% fix nan/inf/neg values to zero
omat(isinf(omat)) = 0;
omat(isnan(omat)) = 0;
omat(omat < 0) = 0;

end
