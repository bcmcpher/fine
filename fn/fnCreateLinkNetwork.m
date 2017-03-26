function [ omat, out, fh ] = fnCreateLinkNetwork(fe, pconn, label)
%% create link network
% Brent McPherson
% 20170212
% 
% Add other measures?
% Cosine Similarity - pdist2(ld1, ld2, 'cosine') of node indices or
% dictionary orientations? both?
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

display('Computing link intersections...');

% for every connections overlap
parfor ii = 1:length(pairs)
    
    % simple index for each pconn combination
    ti1 = pairs(ii, 1);
    ti2 = pairs(ii, 2);
    
    % grab the streamline indices of each link
    li1 = getfield(pconn{ti1}, label, 'indices');
    li2 = getfield(pconn{ti2}, label, 'indices');
    
    % if either connection is empty, fill in 0 and move on
    if isempty(li1) || isempty(li2)
        out{ii}.indices = [];
        out{ii}.dice = 0;
        continue
    end
    
    % create indices of link subtensors
    ln1 = find(fe.life.M.Phi(:,:,li1));
    ln2 = find(fe.life.M.Phi(:,:,li2));
    
    % grab unique voxel indices from each link - 
    % the unique 2nd dimension of the subtensors
    vxln1 = unique(ln1(:, 2));
    vxln2 = unique(ln2(:, 2));
    
    % find the link intersection
    indices = intersect(vxln1, vxln2);
    
    % if the intersection is empty, move on
    if isempty(indices)
        out{ii}.indices = [];
        out{ii}.dice = 0;
        continue
    end
    
    % find the size of the intersection
    num = size(indices, 1);
    
    % combine the voxel indices for denominator of Dice coeff
    den = size(unique([ vxln1; vxln2 ]), 1);
    
    % build Dice coeff values for assignment
    val = (2 * num) / den;
    
    % catch output
    out{ii}.indices = indices;
    out{ii}.dice = val;
    
end

%% assemble matrix

display('Creating link network...');

for ii = 1:length(pairs)
    
    % simple index for each 
    ti1 = pairs(ii, 1);
    ti2 = pairs(ii, 2);
    
    % build matrix from out
    omat(ti1, ti2) = out{ii}.dice;
    omat(ti2, ti1) = out{ii}.dice;
    
end

clear ii ti1 ti2 

% fix nan/inf values to zero
omat(isinf(omat)) = 0;
omat(isnan(omat)) = 0;
%omat(omat < 0) = 0;

%% simple plot    

plt = log10(omat);
%plt = cfdat.agree;

fh = figure();
title('Link Network');
colormap('hot');
imagesc(plt);
axis('square'); axis('equal'); axis('tight');
colorbar;
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);

end
    
%% cosine similarity metric development

% % the data values
% x = pconn{2}.vxend;
% y = pconn{15}.vxend;
% 
% % manually?
% theta = (x * y') / (norm(x, 2) * norm(y, 2));
% theta = acos((x * y') / (norm(x, 2) * norm(y, 2)));
% 
% % with built-in?
% theta = pdist([x; y], 'cosine'); % what do I do with the vector?
