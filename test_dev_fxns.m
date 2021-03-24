%% run the current toolbox to track what functionality is fully working

% load the the data for testing
load_dev_data

% create an initial network output
[ netw, out ] = feCreatePairedConnections(parc, fg, names, 0, 'weights', wght);

% create central tendency measures of FA for each edge
netw = fnAverageEdgeProperty(netw, fg, fa, 'fa');

% create tract profiles of FA for each edge
netw = feTractProfilePairedConnections(netw, fg, fa, 'fa');

% create tract shape data for each edge
netw = fnTractCurvePairedConnections(netw, fg);



