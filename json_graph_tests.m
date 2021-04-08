%% load the json graph data file

% load an example structure
njg = loadjson('/geode2/home/u010/bcmcpher/Carbonate/fine-data/json-graph/network.json');
njn = loadjson('/geode2/home/u010/bcmcpher/Carbonate/fine-data/json-graph/null_network.json');
njs = loadjson('/geode2/home/u010/bcmcpher/Carbonate/fine-data/json-graph/stats_network.json');

% export netw to json
mjg = fnCreateJSONGraph(netw, 'density');

% save and reload my object
savejson('', mjg, '/geode2/home/u010/bcmcpher/Carbonate/fine-data/json-graph/my_network.json');
ljg = loadjson('/geode2/home/u010/bcmcpher/Carbonate/fine-data/json-graph/my_network.json', 'FastArrayParser', 1);

% the basic structure to build
out.graph.label = 'probably the type of edge weight out';
out.graph.directed = 0;

% the network metadata field
out.graph.metadata.unit = 'xyz^3';                 % the space it's it?
out.graph.metadata.desc = 'whatever I want here?'; % depending on what else i guess
out.graph.metadata.extra_name = 'self-loop';       % extra defs?
out.graph.metadata.extra_desc = 'index(x,x) is the diagonal'; % what else can go here?
% anything else I can add?

% node structure names need a string label
out.graph.nodes.('somename').label = '2044';
out.graph.nodes.('somename').metadata.name = 'ParcName';
out.graph.nodes.('somename').metadata.voxel_value = 2044; % basically always label
% check how stats are added

% edges seem to be more similar
out.graph.edges{1}.source = 1; % index of ROI1
out.graph.edges{1}.target = 2; % index of ROI2
out.graph.edges{1}.metatdata.weight = 1; % the edge weight
% stats get appended here?

