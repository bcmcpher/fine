function [ out ] = fnAppNetworkLinkMatrix(config)
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

% set up the link network options
switch config.lnet
    case 'dice'
        ltype = 'dice';
        lcnfg = [];
        lindx = 6;
        ljson = 'volume';
    case 'mi'
        ltype = 'mi';
        lcnfg = config.mname;
        lindx = 7;
        ljson = config.mname;
    otherwise
        error('This is an impossible option.');
end

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
netw = fnAveragePropertyEdges(netw, fg, edata, config.mname, true);

% estimate the link network
netw = fnCreateLinkNetwork(netw, ltype, lcnfg);

% create matrix field for export of all computed edge weights
netw = fnComputeMatrixEdges(netw);

% create connectivity matrices
omat = fnCreateAdjacencyMatrices(netw);

% extract the specific network to write out
nmat = squeeze(omat(:,:,lindx));

% create the link network as a sparse matrix
lmat = fnCreateLinkMatrix(netw, 1);

%% write out multiple predetermined networks

% save the network .mat structure
mkdir(fullfile(pwd, 'raw'));
save(fullfile(pwd, 'raw', 'netw.mat'), 'netw', 'lmat', '-v7.3');

% make nested folders for conmat / network json output for regular network
mkdir(fullfile(pwd, 'conmat-network'));
mkdir(fullfile(pwd, 'conmat-network', 'csv'));
mkdir(fullfile(pwd, 'conmat-links'));
mkdir(fullfile(pwd, 'conmat-links', 'csv'));

% make nested folders for conmat / network json output for link network
mkdir(fullfile(pwd, 'network-json'));

% write the network to the output folder
dlmwrite(fullfile('conmat-network', 'csv', 'volume.csv'), nmat, 'delimiter', ',');

% copy label.json to each output - probably redundant
copyfile('label.json', 'conmat-network');

% write index.json w/ basic descriptors - what is this for?
index = struct('filename', strcat(ltype, '.csv'), 'unit', ltype, 'name', ltype, 'desc', 'FiNE - weighted connectivity');
savejson('', index, fullfile(pwd, 'conmat-network', 'index.json'));

% write jsongraph to disc
fnCreateJsonGraph(netw, ljson, 'network-json/network.json');
gzip('network-json/network.json'); % gzip by defualt
delete('network-json/network.json'); % remove the unzipped file

% write the link network matrix to disc
dlmwrite(fullfile('conmat-links', 'csv', strcat(ltype, '.csv')), full(lmat), 'delimiter', ',');

% build the label.json for the link network

% pull the number of nodes to rebuild possible edge list
nnodes = size(olabel, 1);

% build the upper diagonal node indices
nmsk = triu(ones(nnodes), 1);
[ xnidx, ynidx ] = find(nmsk > 0);

% build the possible edges that could exist
eindx = sortrows([ xnidx, ynidx ]);
nlink = size(eindx, 1);
clear nnodes nmsk xnidx ynidx

% build a edge by edge set of labels
llabel = cell(nlink, 1);
for link = 1:nlink
    llabel{link}.name = [ olabel{eindx(link,1)}.name '-' olabel{eindx(link,2)}.name ];
    llabel{link}.label = nan;
    llabel{link}.voxel_value = [ olabel{eindx(link,1)}.voxel_value, olabel{eindx(link,2)}.voxel_value ];
end

% write edge label.json to to disk
savejson('', llabel, 'conmat-links/label.json');

% write index.json w/ basic descriptors - what is this for?
index = struct('filename', strcat(ltype, '.csv'), 'unit', ltype, 'name', ltype, 'desc', 'FiNE - link network');
savejson('', index, fullfile(pwd, 'conmat-links', 'index.json'));

% there's the 2 ways to deal with exporting link networks to json-graph, 
% and it's not obvious which one is worth the effort yet. There's no use
% anyway, so wait.

out = 'done.';

end
