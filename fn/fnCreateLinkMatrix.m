function [ omat ] = fnCreateLinkMatrix(netw, prop)
%[ omat ] = fnCreateLinkMatrix(netw, prop);
%   Create a sparse matrix containing the non-zero elements estimated from
%   fnCreateLinkNetwork(). Creates only a single property.
%
%   Optionally plot?
%

%% parse input

if ~isfield(netw, 'links')
    error('No link network is computed.');
end

% if dtype isn't passed, set as empty b/c it isn't used
if(~exist('prop', 'var') || isempty(prop))
    prop = 1;
end

% add 2 for numeric index to look at first sparse value stored
index = prop + 2;

% if the index is in excess of the dimensions, throw a warning and 
% reset to the first valid index.
if ((size(netw.links.values, 2)) < index)
    warning('The requested numeric index exceeds the number of stored dimensions. Resetting to the first value.');
    index = 3;
end

%% pull data

% compute the size of the link network
nnet = triu(ones(size(netw.nodes, 1)), 1); % pull the upper diagonal
mdim = sum(nnet(:) > 0); % pull the possible number of edges

% pull the data from the netw links
data = netw.links.values;

%% stack the data to create a symmetric sparse matrix

disp([ 'Creating link network graph with: ' netw.links.labels{index} ]);

% pull the indices, stack for symmetry (i->j == j->i)
idx = [ data(:,1), data(:,2); data(:,2), data(:,1) ];

% pull the property, stack for symmetry
val = [ data(:, index); data(:, index) ];

% create a symmetric sparse matrix from stored data
omat = sparse(idx(:,1), idx(:,2), val, mdim, mdim);

end
