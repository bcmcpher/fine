function [ omat, out ] = fnCreateLinkNetwork(pconn, label)
%% create link network
% Brent McPherson
% 20170212
%
% Need to figure out tensor indice to ROI coords mapping
%
% Add other measures?
% Cosine Similarity - pdist(ld, 'cosine'); or pdist2(ld1, ld2, 'cosine'); 
% where ld is matrix of node indices, counts, and/or dictionary orientations.
% some of these? all of these?
%
% https://en.wikipedia.org/wiki/Diversity_index#Simpson_index
% Simpson's Index: probability that 2 random voxels are part of intersection
% Richness: how many links does an intersection contain
% Shannon Diversity Index: global measure of intersections between links 
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
    
    % grab the unique voxels of each edge
    li1 = getfield(pconn{ti1}, label, 'pvoxels');
    li2 = getfield(pconn{ti2}, label, 'pvoxels');
    
    % if either connection is empty, fill in 0 and move on
    if isempty(li1) || isempty(li2)
        out{ii}.indices = [];
        out{ii}.dice = 0;
        continue
    end
    
    % find the link intersection
    %indices = intersect(li1, li2, 'rows');
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
    %den = size(unique([ li1; li2 ]), 1); % values greater than 1
    den = size(li1, 1) + size(li2, 1); % values less than 1
    
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
