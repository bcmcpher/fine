function [ omat, outLabels ] = feCreateAdjacencyMatrices(pconn, label)
%feCreateAdjacencyMatrices create adjacency matrices from a pconn list and
% returns the labels for each matrix along the third dimension.
%   

% build matrices from pconn

display('Building Adjacency Matrices...');

% pull .matrix component from cell to determine how many matrices to make
tmp = getfield(pconn{1}, label);
tmp = tmp.matrix;

% find the and label fields
connLabels = fieldnames(tmp);
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
        val =  getfield(pconn{ii}, label);
        val = val.matrix;
        val = getfield(val, connLabels{jj});
        
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
