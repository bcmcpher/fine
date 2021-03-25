%% run the current toolbox to track what functionality is fully working

% load the the data for testing
load_dev_data

% create an initial network output
[ netw, out ] = fnCreateEdges(lobe, fg, lname, 0, 'weights', wght);
%[ netw, out ] = fnCreateEdges(node, fg, nname, 0, 'weights', wght);

% create central tendency measures of FA for each edge
netw = fnAverageEdgeProperty(netw, fg, fa, 'fa');

% create tract profiles of FA for each edge
netw = feTractProfilePairedConnections(netw, fg, fa, 'fa');

% create tract shape data for each edge
netw = fnTractCurvePairedConnections(netw, fg);

% run virtual lesion
netw = feVirtualLesionPairedConnections(netw, M, wght, dsig, nTheta, S0);

% create matrix field for export of all computed edge weights
netw = fnComputeMatrixField(netw);

% create connectivity matrices
[ omat, olab ] = feCreateAdjacencyMatrices(netw);

% check the count matrix
figure; imagesc(omat(:,:,1)); title(olab{1}); colorbar;
