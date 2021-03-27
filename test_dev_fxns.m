%% run the current toolbox to track what functionality is fully working

% load the the data for testing
load_dev_data

% create an initial network output
[ netw, out ] = fnCreateEdges(lobe, fg, lname, 0, 'weights', wght);
%[ netw, out ] = fnCreateEdges(node, fg, nname, 0, 'weights', wght);

% create central tendency measures of FA for each edge
netw = fnAveragePropertyEdges(netw, fg, fa, 'fa');

% create tract profiles of FA for each edge
netw = fnTractProfileEdges(netw, fg, fa, 'fa');

% create tract shape data for each edge
netw = fnTractCurveEdges(netw, fg);

% run virtual lesion
netw = fnVirtualLesionEdges(netw, M, wght, dsig, nTheta, S0);

% create matrix field for export of all computed edge weights
netw = fnComputeMatrixEdges(netw);

% create connectivity matrices
[ omat, olab ] = fnCreateAdjacencyMatrices(netw);

% check the count matrix
figure; imagesc(omat(:,:,1)); title(olab{1}); colorbar;

%% scratch

% dice coefficient between orig / acpc tree assignment
dcoeff = zeros(length(net1.edges), 1);
for ii = 1:length(net1.edges)
    dcoeff(ii) = 2*length(intersect(net1.edges{ii}.fibers.indices, net2.edges{ii}.fibers.indices)) / (length(net1.edges{ii}.fibers.indices) + length(net2.edges{ii}.fibers.indices));
end
% rather different, but probably more robust?

% if anatomy looks good, keep it.

