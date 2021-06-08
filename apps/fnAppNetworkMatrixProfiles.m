function [ out ] = fnAppNetworkMatrixProfiles(config)
%[ out ] = fnAppNetworkAverageEdgeMatrix(config);
%   The version of this tool to run on brainlife. Needs to be split up into
%   many apps b/c the platform favors lots of little fxns over one flexible
%   one. 
%
%   This loads a 'config.json' file with sufficiently valid inputs and
%   writes a brainlife.io conmat data type and a jsongraph (json.gz) data
%   type for the requested average property of the edge as the weights 
%   estimated for a structural network.
%
% Brent McPherson (c) 2021, Indiana University
%

%% parse the config.json inputs

% load the config.json
config = loadjson(config);

% load .tck streamlines without downsampling
fg = dtiImportFibersMrtrix(config.track, 0.5);

% load parcellation data
parc = niftiRead(config.parc);

% load the passed labels
jlabel = loadjson(config.label)';

% pull logical for cleaning edges or not
clean_edges = config.clean;

% pull logical for getting edge shapes
edge_shapes = config.shapes;

%% fix label.json values

% I don't try to parse the original labels - if you need to refer to the original 
% labels from whatever reference uses them then use a tool that doesn't relabel the 
% voxels values so they no longer correspond to the labels for no clear reason.
%label = cellfun(@(x) str2num(x.label), jlabel);

% pull names and remove irrelevant stems / illegal characters for fieldnames
names = cellfun(@(x) { x.name x.voxel_value }, jlabel, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '.label', ''), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_LH_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_RH_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '+', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '.', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_L_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_R_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_ROI', ''), x{2} }, names, 'UniformOutput', false);
% there will be more bad naming conventions to fix here I'm sure...

% pull the values in voxels present in the volume
voxel = cellfun(@(x) x.voxel_value, jlabel);

% pull the unique non-zero values present in the voxels
uparc = unique(parc.data(:));
uparc = uparc(uparc > 0);

% get the values present in the parc.nii.gz that are missing in the labels.json
lmiss = setdiff(uparc, voxel);

% deal with incomplete maTT labels / still be robust to wrong labels passed
if ~isempty(lmiss) % is there is any mismatch in parc/json
    
    % look for maTT bug
    if size(lmiss, 1) == 14
        disp('The passed label.json is missing 14 observations. This is a common issue with maTT.');
        disp('Appending ''names'' with what the missing labels are assumed to be.');
        
        % pull the labels - assume 0-max cortical
        mlab = max(cellfun(@(x) x{2}, names));
        
        % build the subcortical names
        lhsub = {'lh_Thalamus'; 'lh_Caudate'; 'lh_Putamen'; 'lh_Pallidum'; 'lh_Hippocampus'; 'lh_Amygdala'; 'lh_Accumbens'};
        rhsub = strrep(lhsub, 'lh_', 'rh_');
        wbsub = [ lhsub; rhsub ];
        
        % build the subcortical name/values
        subnm = cellfun(@(x,y) {x, y}, wbsub, num2cell([mlab+1:mlab+14]'), 'UniformOutput', false);
        
        % append the name/value combination
        names = [ names; subnm ];

        % create a fixed labels.json because the conmat type needs it
        olabel = cell(size(uparc, 1), 1);
        for ii = 1:size(uparc, 1)
            olabel{ii}.name = names{ii}{1}; % this assumes they're sorted, but they are in maTT
            olabel{ii}.label = '???';    % b/c the relabeling enforces it.
            olabel{ii}.voxel_value = names{ii}{2};
        end
        
        % write fixed label.json to to disk
        savejson('', olabel, 'label.json');

    else % otherwise it's just wrong
        error('The passed label.json has does not match the values stored in parc.nii.gz');
    end
    
else
    
    % otherwise the provided labels actually match the parcellation like
    % they are supposed to - so just copy it over
    copyfile(config.label, '.');
    
end

clear jlabel lmiss ii lhsub rhsub

%% determine if the requested edge measure is actually passed as input

% check requested name across all possible inputs for path to file to load
switch config.mname
    % tensor
    case 'fa'
        efile = config.fa;
    case 'md'
        efile = config.md;
    case 'rd'
        efile = config.rd;
    case 'ad'
        efile = config.ad;
    case 'cl'
        efile = config.cl;
    case 'cp'
        efile = config.cp;
    case 'cs'
        efile = config.cs;
    % kurtosis
    case 'ga'
        efile = config.ga;
    case 'mk'
        efile = config.mk;
    case 'ak'
        efile = config.ak;
    case 'rk' 
        efile = config.rk;        
    % noddi
    case 'ndi'
        efile = config.ndi;
    case 'isovf'
        efile = config.isovf;
    case 'odi'
        efile = config.odi;
    % myelin map
    case 'map'
        efile = config.map;
    otherwise
        error('The requested file is not configured in this input: %s ', config.mname);
end
% technically, any volume in register to the input data can be sampled here

% try to load the requested file
try
    edata = niftiRead(efile);
catch
    error('The requested output - %s - is not stored for this input.', config.mname);
end

%% actually build the networks after fixing labels

% create the network object
netw = fnCreateEdges(parc, fg, names, config.maxDist, config.minStrm);

% optionally clean the edges
if clean_edges
    netw = fnCleanEdges(netw, fg, config.minLength, config.maxVars, config.maxLengthStd, config.numNodes, config.maxIter, config.minStrm);
end

% estimate the average edge property
netw = fnTractProfileEdges(netw, fg, edata, config.mname, config.numNodes, config.minStrm);

% % optionally get the edge shape data
% if edge_shapes
%     netw = fnTractCurveEdges(netw, fg, config.nnodes, config.minNum);
% end

%% write out the profiles

% pull the node labels
nodes = cellfun(@(x) x.name, netw.nodes, 'UniformOutput', false);
edgec = size(netw.edges, 1);

% make the profile output
mkdir('./tractmeasures');

% open the output file
fileID = fopen('tractmeasures/tractmeasures.csv','w');

% put the first label
fprintf(fileID, '%s,', 'structureID');

% print the column headers
for node = 1:config.numNodes
    if node == config.numNodes
        fprintf(fileID, '%s', sprintf('node%03d', node));
    else
        fprintf(fileID, '%s,', sprintf('node%03d', node));
    end
end

% line break
fprintf(fileID, '\n');

% write out tractmeasures.csv file
for edge = 1:edgec
    
    % write a label to file
    structureID = strcat(nodes{netw.parc.pairs(edge,1)}, '-', nodes{netw.parc.pairs(edge,2)});
    fprintf(fileID,'%s,', structureID{1});
    
    % pull the profile
    tprof = netw.edges{edge}.profile.(config.mname);
    
    % for every node in a profile
    for node = 1:config.numNodes
        
        % store the nodes
        if node == config.numNodes
            fprintf(fileID,'%.4f\n', tprof(node));
        else
            fprintf(fileID,'%.4f,', tprof(node));
        end
        
    end
    
end

% close the output file
fclose(fileID);

out = 'done.';

end
