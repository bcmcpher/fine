function [ omat, outLabels ] = feCreateAdjacencyMatrices(pconn, label)
%feCreateAdjacencyMatrices creates adjacency matrices from a pconn array and
% returns the labels for each matrix along the third dimension. The number
% of matrices returned depends on the number of fields computed for pconn.
%   
% INPUTS:
%     pconn - is the paired connections object to create adjacency matrices from.
%
%     label - string indicating the fiber groups for which to create virtual lesions
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%             Additionally, this can be run after cleaning, resulting in
%             valid calls of 'all_clean' and 'nzw_clean', respectively.
%
% OUTPUTS:
%     omat - 3d array containing (nodes x nodes x edge_type) of processed networks
%     olab - cell array of labels of the 'edge_type' along the 3rd dimension of omat
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
% Brent McPherson (c), 2017 - Indiana University
%

% build matrices from pconn

display('Building Adjacency Matrices...');

% pull .matrix component from cell to determine how many matrices to make
% this assumes all the fields are the same - reasonably true
tmp = pconn{1}.(label);
tmp = tmp.matrix;

% find the and label fields
connLabels = fieldnames(tmp);
outLabels = cell(length(connLabels), 1);
for ii = 1:length(connLabels)
    outLabels{ii} = [ label '_' connLabels{ii} ];
end

% find output size
nmat = size(connLabels, 1);

% find the number of unique labels
uniquelabels = zeros(size(pconn, 1), 1);
for ii = 1:size(pconn, 1)
    uniquelabels(ii, 1) = pconn{ii}.roi1;
end
nlabs = size(unique(uniquelabels), 1) + 1;

% initialize output
omat = zeros(nlabs, nlabs, nmat);

% re-create pairs based on number of nodes
pairs = nchoosek(1:nlabs, 2);

clear tmp

% dummy check until I can think better
if size(pairs, 1) ~= size(pconn, 1)
    error('The paired connections do not match the expected number of ROI pairings.');
end

% for every paired connection
for ii = 1:length(pconn)
    for jj = 1:nmat
        
        % catch the specific matrix value
        val = pconn{ii}.(label);
        val = val.matrix;
        val = val.(connLabels{jj});
        
        % assign labels to output matrix
        omat(pairs(ii, 1), pairs(ii, 2), jj) = val;
        omat(pairs(ii, 2), pairs(ii, 1), jj) = val;
        
    end
end

% fix impossible values
omat(isnan(omat)) = 0;
omat(isinf(omat)) = 0;
omat(omat < 0) = 0;

end
