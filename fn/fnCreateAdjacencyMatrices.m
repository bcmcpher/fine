function [ omat, olab ] = fnCreateAdjacencyMatrices(netw, srt)
%feCreateAdjacencyMatrices() creates adjacency matrices from a pconn array and
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
%     srt   - (optional) apply a specific order to the cortical nodes to group the
%             connections differently. By default, the nodes are sorted by
%             the number assigned in the cortical label volume.
%
% OUTPUTS:
%     omat - 3d array containing (nodes x nodes x edge_type) of processed networks
%
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

% parse optional arguments
if(~exist('srt', 'var') || isempty(srt))
    srt = [];
end

% build matrices from pconn

disp('Building Adjacency Matrices...');

% pull .matrix component from cell to determine how many matrices to make
% this assumes all the fields are the same - reasonably true
tmp = netw.edges{1}.matrix;

% find the and label fields
olab = fieldnames(tmp);

% find output size
nmat = size(olab, 1);

% pull the number of labels
nlabs = size(netw.nodes, 1);

% initialize output
omat = zeros(nlabs, nlabs, nmat);

% re-create pairs based on number of nodes
pairs = netw.parc.pairs;

clear tmp

% for every paired connection
for ii = 1:length(netw.edges)
    for jj = 1:nmat
        
        % catch the specific matrix value
        val = netw.edges{ii}.matrix.(olab{jj});
        
        % assign labels to output matrix
        omat(pairs(ii, 1), pairs(ii, 2), jj) = val;
        omat(pairs(ii, 2), pairs(ii, 1), jj) = val;
        
    end
end

% fix impossible values - kick a warning?
omat(isnan(omat)) = 0;
omat(isinf(omat)) = 0;
%omat(omat < 0) = 0;

% sort the cortical nodes given a vector
if ~isempty(srt)
    disp('The output matrices are sorted by the user provided node order.');
    omat = omat(srt, srt, :);
end

end
