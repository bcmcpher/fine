function [ omat, olab, out ] = fnCreateLinkNetwork(pconn, label, feroi, img)
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
%                     'zwr' for all zero weightef fibers removed by LiFE
%             Additionally, this can be run after cleaning, resulting in
%             valid calls of 'all_clean' and 'nzw_clean', respectively.
%
% OUTPUTS:
%     omat - a matrix that is (length(pconn) x length(pconn)) in size
%            showing the Dice coefficient of each intersecting volume.
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
% vx            = [ 1.25 1.25 1.25 ];
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % find the voxels that each edge occupies
% pconn = fnFindPathVoxels(pconn, 'nzw', Phi, vx);
%
% % create link network based on overlap of non-zero weighted fibers
% [ omat, out ] = fnCreateLinkNetwork(pconn, 'nzw');
%
% Brent McPherson (c), 2017 - Indiana University
%

%% generate link network

% error if volume is not precomputed
if (~isfield(pconn{1}.(label), 'volume'))
    error('Edge volumes for this label must be precomputed with fnFindPathVoxels().');
end    

display('Preallocating output structure...');

% number of edges; link nodes
lnodes = size(pconn, 1);

% create indices for every pair of link connections
pairs = nchoosek(1:lnodes, 2);

% build empty cells to parallelize computation
out = cell(length(pairs), 1);

% build empty matrices
omat = zeros(lnodes, lnodes, 3);

display('Computing unique edge intersections...');

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
    li1 = pconn{ti1}.(label).volume.pvoxels;
    li2 = pconn{ti2}.(label).volume.pvoxels;
    
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
    
    % BEGIN PULLING VALUES FOR JOINT ENTROPY / MI HERE
    % THESE ARE ONLY CONNECTIONS WITH SHARED SPACES

    %keyboard;
    
    % grab image size and ROI image space coordinates
    imc1 = feroi(li1, :);
    imc2 = feroi(li2, :);

    % preallocate image values
    dims = max([size(imc1, 1), size(imc2, 1)]);
    im1 = zeros(dims, 1);
    im2 = zeros(dims, 1);
    
    % because I can't figure out easy indexing, loop to get each voxel
    % intensity for both images
    for jj = 1:size(imc1, 1)
        im1(jj) = img.data(imc1(jj, 1), imc1(jj, 2), imc1(jj, 3));
    end
    
    for jj = 1:size(imc2, 1)
        im2(jj) = img.data(imc2(jj, 1), imc2(jj, 2), imc2(jj, 3));
    end
        
    % compute bins of joint histogram
    [ ~, ~, indrow ] = unique(im1(:));
    [ ~, ~, indcol ] = unique(im2(:));
    
    % compute joint entropy
    jointHistogram = accumarray([indrow indcol], 1);
    jointProb = jointHistogram / numel(indrow);
    indNoZero = jointHistogram ~= 0;
    jointProb1DNoZero = jointProb(indNoZero);
    jointEntropy = -sum(jointProb1DNoZero.*log2(jointProb1DNoZero));
    
    % compute individual histogram summaries
    histogramImage1 = sum(jointHistogram, 1);
    histogramImage2 = sum(jointHistogram, 2);
    
    % find non-zero elements for first image's histogram
    indNoZero = histogramImage1 ~= 0;
    
    % extract them out and get the probabilities
    prob1NoZero = histogramImage1(indNoZero) / numel(histogramImage1);
    
    % compute the entropy
    entropy1 = -sum(prob1NoZero.*log2(prob1NoZero));
    
    % repeat for the second image
    indNoZero = histogramImage2 ~= 0;
    prob2NoZero = histogramImage2(indNoZero) / numel(histogramImage2);
    entropy2 = -sum(prob2NoZero.*log2(prob2NoZero));
    
    % now compute mutual information
    mutualInformation = entropy1 + entropy2 - jointEntropy;
    
    % END OF IMAGE INTENSITY COMPARISONS
    
    % find the size of the intersection
    num = size(indices, 1);
    
    % combine the voxel indices for denominator of Dice coeff
    den = size(li1, 1) + size(li2, 1);
    
    % build Dice coeff values for assignment
    val = (2 * num) / den;
    
    % catch output
    out{ii}.indices = indices;
    out{ii}.dice = val;
    out{ii}.mutualInformation = mutualInformation;
    out{ii}.jointEntropy = jointEntropy;
    out{ii}.histogramImage1 = histogramImage1;
    out{ii}.histogramImage2 = histogramImage2;
   
end
time = toc;

display(['Created link network from ' num2str(length(pairs)) ' unique edge combinations in ' num2str(round(time)/60) ' minutes.']);

%% assemble matrix

display('Creating link network matrix...');

for ii = 1:length(pairs)
    
    % simple index for each 
    ti1 = pairs(ii, 1);
    ti2 = pairs(ii, 2);
    
    % build matrix of dice coeffs
    omat(ti1, ti2, 1) = out{ii}.dice;
    omat(ti2, ti1, 1) = out{ii}.dice;
    
    % build matrix of dice coeffs
    omat(ti1, ti2, 2) = out{ii}.jointEntropy;
    omat(ti2, ti1, 2) = out{ii}.jointEntropy;
    
    % build matrix of dice coeffs
    omat(ti1, ti2, 3) = out{ii}.mutualInformation;
    omat(ti2, ti1, 3) = out{ii}.mutualInformation;
    
end

clear ii ti1 ti2 

% fix nan/inf/neg values to zero
omat(isinf(omat)) = 0;
omat(isnan(omat)) = 0;
omat(omat < 0) = 0;

olab = {'dice', 'jointEntropy', 'mutualInformation'};

end
