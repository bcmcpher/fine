function [ jsg ] = fnCreateJSONGraph(netw, edgew)
%[ jsg ] = fnCreateJSONGraph(netw, edgew);
%   Write out a .json graph object of the network structure.
%
% INPUTS:
%     netw  - the network object with matrix fields estimated
%     edgew - the edge weight estimated in the network to store
%
% OUTPUTS:
%     jsg - the json graph structure to save with 
%
% Brent McPherson (c) 2021, Indiana University
%

%% parse the arguments, make sure the requested field exists

% set a default maximum distance in mm an endpoint can be from a node.
if(~exist('edgew', 'var') || isempty(edgew))
    edgew = 'density';
end

% if the requested field doesn't exist
if ~isfield(netw.edges{1}, 'matrix')
    error('Matrix field must be computed with fnComputeMatrixEdges();');
else
    if ~isfield(netw.edges{1}.matrix, edgew)
        warning('The requested edge weight to write out does not exist. Returning ''density'' by default.');
        edgew = 'density';
    end
end

%% initialize the graph structure / sizes

% pull the dimensions of nodes / edges
nnodes = length(netw.nodes);
nedges = length(netw.edges);

% the basic structure to build
jsg.graph.label = edgew;
jsg.graph.directed = 0;

% the network metadata field
jsg.graph.metadata.unit = edgew; % the unit of values in the edges
jsg.graph.metadata.desc = [ 'FiNE - ' edgew ' weighted connectivity' ]; 
% any naming convention?

%jsg.graph.metadata.extra_name = 'self-loop';   % extra defs?
%jsg.graph.metadata.extra_desc = 'index(x,x) is the diagonal'; % what else can go here?
% anything else I can add?

% how would stats be added?

%% initialize and store the node data

% preallocate struct array of nodes w/ names?
nmetad = struct('name', [], 'voxel_value', [], 'volume', [], 'center', []);
tnode = struct('label', [], 'metadata', nmetad);

% append node information as struct array
for node = 1:nnodes
    
    % sanitize stored labels for valid output fields
    label = strrep(netw.nodes{node}.name{1}, '-', '_');
    
    % minimum data for node
    nodes.(label) = tnode;
    nodes.(label).label = netw.nodes{node}.label;
    nodes.(label).metadata.name = netw.nodes{node}.name;
    nodes.(label).metadata.voxel_value = netw.nodes{node}.label;
    % any additional fields?
    
    % my additions of node data
    nodes.(label).metadata.volume = netw.nodes{node}.volume;
    nodes.(label).metadata.center = netw.nodes{node}.center.acpc;

    % check how stats would be added

end

% add graph to jsg
jsg.graph.nodes = nodes;

%% store the edge data

% append cell information as a cell array
for edge = 1:nedges
    
    % store the minimum data for each edge
    jsg.graph.edges{edge}.source = netw.parc.pairs(edge, 1); 
    jsg.graph.edges{edge}.target = netw.parc.pairs(edge, 2);
    jsg.graph.edges{edge}.metatdata.weight = netw.edges{edge}.matrix.(edgew);
    % only 1 edge weight - they can't handle multiple weights
    
    % other info to store here?
    
    % stats get appended here?

end

end
