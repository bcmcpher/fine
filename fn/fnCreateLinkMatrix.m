function [ omat ] = fnCreateLinkMatrix(netw, prop)
%[ omat ] = fnCreateLinkMatrix(netw, prop);
%   Create a sparse matrix containing the non-zero elements estimated from
%   fnCreateLinkNetwork(). Creates only a single property.
%
%   Optionally plot?
%

%% pull basic data

% pull the size of the link network
mdim = size(netw.edges, 1);

% pull the data from the network
data = netw.links.values;

%% check that prop index is a valid call

if prop < 0
    error('You must request a positive value. There are no negative indices.');
end

if prop < 2
    warning('The first 2 properties store the edge indices. Setting ''prop'' to first valid index (3).');
    prop = 3;
end

if prop > size(data, 2)
    error('The requested property only has %d values stored - ''prop'' is too large.', size(data, 2));
end

%% stack the data to create a sparse matrix

disp([ 'Creating link network graph with: ' netw.links.labels{prop} ]);

% pull the indices, stack for symmetry
idx = [ data(:,1), data(:,2); data(:,2), data(:,1) ];

% pull the property, stack for symmetry
val = [ data(:, prop); data(:, prop) ];

% create a symmetric sparse matrix from stored data
omat = sparse(idx(:,1), idx(:,2), val, mdim, mdim);

end
