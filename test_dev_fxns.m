%% run the current toolbox to track what functionality is fully working

% load the the data for testing
load_dev_data

% create an initial network output
%[ netw, out ] = fnCreateEdges(lobe, fg, lname, 2, 0, 'weights', wght);
[ netw, out ] = fnCreateEdges(node, fg, nname, 4, 5, 'weights', wght);

% create central tendency measures of FA for each edge
netw = fnAveragePropertyEdges(netw, fg, fa, 'fa', true, true, true);

% create tract profiles of FA for each edge
netw = fnTractProfileEdges(netw, fg, fa, 'fa');
netw = fnTractProfileEdges(netw, fg, dt6, 'dt6');

% store the tract profile tensor and average it into modules
tpmat = fnTractProfileTensor(netw, 'fa'); 
[ mnten, seten ] = fnTractProfileModules(tpmat, modules);

% create tract shape data for each edge
netw = fnTractCurveEdges(netw, fg);

% run virtual lesion
netw = fnVirtualLesionEdges(netw, M, wght, dsig, nTheta, S0);

% virtual lesion network modules
[ vmat, ~, indx, pcmodl ] = fnModuleVirtualLesion(netw, modules, modName, M, wght, dsig, nTheta, S0);

% create matrix field for export of all computed edge weights
netw = fnComputeMatrixEdges(netw, 'weights', [], true);

% create connectivity matrices
[ omat, olab ] = fnCreateAdjacencyMatrices(netw);

% test jsongraph out
fnCreateJsonGraph(netw, 'density', 'test.json');

% estimate modules of the network
[ Mden, Mden_in, Mden_bw, Mden_in_lo, Mden_bw_hi ] = fnModuleDensity(squeeze(omat(:,:,1)), modules, 'mean');

% estimate the link network
netw = fnCreateLinkNetwork(netw, 'dice');

% check the count matrix
figure; imagesc(omat(:,:,1)); title(olab{1}); colorbar;

%% scratch

tfg = fg.fibers(netw.edges{5}.fibers.indices);
figure; hold on;
for ii = 1:size(tfg, 1)
    plot3(tfg{ii}(1,:), tfg{ii}(2,:), tfg{ii}(3,:), 'red');
end
hold off; axis image
clear tfg ii

% dice coefficient between orig / acpc tree assignment
dcoeff = zeros(length(net1.edges), 1);
for ii = 1:length(net1.edges)
    dcoeff(ii) = 2*length(intersect(net1.edges{ii}.fibers.indices, net2.edges{ii}.fibers.indices)) / (length(net1.edges{ii}.fibers.indices) + length(net2.edges{ii}.fibers.indices));
end
% rather different, but probably more robust?

% if anatomy looks good, keep it.
tfg = fg.fibers(nett.edges{2}.fibers.indices);
figure; hold on
for ii = 1:size(tfg, 1)
    plot3(tfg{ii}(1,:), tfg{ii}(2,:), tfg{ii}(3,:), 'red');
end
hold off; axis image
clear tfg ii
