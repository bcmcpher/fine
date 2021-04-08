function [ out ] = fnAppNetworkCountMatrix(config)
%[ out ] = fnAppNetworkMatrix(config);
%   The version of this tool to run on brainlife. Needs to be split up into
%   many apps b/c the platform favors lots of little fxns over one flexible
%   one.
%
% Brent McPherson (c) 2021, Indiana University
%

% load the config.json
config = loadjson(config);

% load .tck streamlines without downsampling
fg = dtiImportFibersMrtrix(config.tracks, 0.5);

% load parcellation data
node = niftiRead(config.parc);

% reconcile the labels stored in whatever is passed
% deal w/ maTT not passing the full set of labels actually stored in the volume.
jindex = readjson(config.index);
jnames = readjson(config.names);
% deal with whatever this nonsense is - correcly parse and create valid names

% create the network object
netw = fnCreateEdges(node, fg, names, config.search, config.min);

% create matrix field for export of all computed edge weights
netw = fnComputeMatrixEdges(netw);

% create connectivity matrices
[ omat, olab ] = fnCreateAdjacencyMatrices(netw);

% use selected edge weight to save from the allowed options
switch config.edgw
    case 'count'
        oindx = 1;
    case 'density'
        oindx = 2;
    case 'length'
        oindx = 3;
    case 'denlen'
        oindx = 4;
    otherwise
        error('An impossible option was passed through the config.');
end

% write to conmat datatype layout
mkdir csv
dlmwrite([ './csv/' olab{oindx} '.csv' ], omat(:,:,oindx), 'delimiter', ',');
% copy index.json/names.json to wherever they need to be

% create json graph object
jsg = fnCreateJSONGraph(netw, config.edgew);

% write to disk
% does this go in a folder or not?
savejson('', jsg, './network.json');
gzip('./network.json');

out = 'done.';

end

