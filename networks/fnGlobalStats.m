function [ out ] = fnGlobalStats(pconn, label)
%fnGlobalStats returns some basic summary numbers for a pconn network object
%   
%

% total fibers found
tfib = zeros(size(pconn, 1), 1);

% for every entry
for ii = 1:length(pconn)
    
    % pull the desired label
    tmp = getfield(pconn{ii}, label);
    
    % NEED TO COMPUTE VALUES ON EACH MATRIX VALUE?
    
    % pull the streamline count
    tfib(ii) = size(tmp.indices, 1);
    
end

% graph density estimate
nzvals = tfib > 0;

% pull the number of nodes
uniquelabels = zeros(size(pconn, 1), 1);
for ii = 1:size(pconn, 1)
    uniquelabels(ii, 1) = pconn{ii}.roi1;
end
nlabs = size(unique(uniquelabels), 1) + 1;

% compute graph density
out.graph_density = 2 * sum(nzvals) / nlabs**2;
out.count = sum(tfib);


end

