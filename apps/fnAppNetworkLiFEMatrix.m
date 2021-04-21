function [ out ] = fnAppNetworkLiFEMatrix(config)
%[ out ] = fnAppNetworkAverageEdgeMatrix(config);
%   The version of this tool to run on brainlife. Needs to be split up into
%   many apps b/c the platform favors lots of little fxns over one flexible
%   one. 
%
%   This loads a 'config.json' file with sufficiently valid inputs and
%   writes a brainlife.io conmat data type and a jsongraph (json.gz) data
%   type for the LiFE virtual lesion of the edge as the weights estimated 
%   for a structural network.
%
% Brent McPherson (c) 2021, Indiana University
%

%% parse the config.json inputs

% load the config.json
config = loadjson(config);

% load .tck streamlines without downsampling
disp('Loading saved LiFE model...');
load(config.fe);

% load fibers from fg
fg = feGet(fe, 'fibers acpc');

% pull the fiber weights
weights = feGet(fe, 'fiber weights');

% pull the model
M = feGet(fe, 'model');

% subset to only positively weighted fibers
nonzero = weights > 0;
nonzidx = find(nonzero); % because the sparse tensor can't be logically indexed
fg.fibers = fg.fibers(nonzero);
wght = weights(nonzero);
M.Phi = M.Phi(:,:,nonzidx);

% pull the number of directions
nTheta = feGet(fe, 'nbvals');

% pull the demeaned diffusion signal
dsig = feGet(fe, 'dsigdemeaned by voxel');

% pull the mean diffusion signal
S0 = feGet(fe, 'b0signalimage');

% pull the rmse volume for an average edge weight
rmse = fnRmseVolume(fe);
rmse.fname = 'rmse.nii.gz';

% clear fe for memory usage
clear fe nonzero nonzidx

% load parcellation data
disp('Loading parcellation...');
parc = niftiRead(config.parc);

% load the passed labels
jlabel = loadjson(config.label)';

%% fix label.json values

% I don't try to parse the original labels - if you need to refer to the original 
% labels from whatever reference uses them then use a tool that doesn't relabel the 
% voxels values so they no longer correspond to the labels for no clear reason.
%label = cellfun(@(x) str2num(x.label), jlabel);

% pull names and remove irrelevant stems / illegal characters for fieldnames
names = cellfun(@(x) x.name, jlabel, 'UniformOutput', false);
names = strrep(names, '.label', '');
names = strrep(names, '_LH_', '_');
names = strrep(names, '_RH_', '_');
names = strrep(names, '+', '_');
names = strrep(names, '.', '_');
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
        lhsub = {'lh_Thalamus'; 'lh_Caudate'; 'lh_Putamen'; 'lh_Pallidum'; 'lh_Hippocampus'; 'lh_Amygdala'; 'lh_Accumbens'};
        rhsub = strrep(lhsub, 'lh_', 'rh_');
        names = [ names; lhsub; rhsub ];
        
        % create a fixed labels.json because the conmat type needs it
        olabel = cell(size(uparc, 1), 1);
        for ii = 1:size(uparc, 1)
            olabel{ii}.name = names{ii}; % this assumes they're sorted, but they are in maTT
            olabel{ii}.label = '???';    % b/c the relabeling enforces it.
            olabel{ii}.voxel_value = ii;
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

%% actually build the networks after fixing labels

% create the network object
netw = fnCreateEdges(parc, fg, names, config.maxDist, config.minStrm, 'weights', wght);

% estimate the average RMSE within each edge
netw = fnAveragePropertyEdges(netw, fg, rmse, 'rmse', false);

% run virtual lesion
netw = fnVirtualLesionEdges(netw, M, wght, dsig, nTheta, S0);

% create matrix field for export of all computed edge weights
netw = fnComputeMatrixEdges(netw);

% create connectivity matrices
[ omat, olab ] = fnCreateAdjacencyMatrices(netw);

%% write out multiple predetermined networks

% use the default names in olab b/c they're always the same outs
edgei = [ 1 2 4 5 6 7 8 9 12 14 16 ];
edgew = olab(edgei);
edgec = 11;

% write out each network
for edge = 1:edgec
    
    % create folder names for outputs b/c brainlife can't parse paths to outputs
    conout = [ 'conmat-', edgew{edge} ];
    netout = [ 'network-', edgew{edge} ];

    % make a folder in conmat for each modality
    mkdir(fullfile(pwd, conout));
    mkdir(fullfile(pwd, conout, 'csv'));
    
    % write the csv of the adjacency matrix (conmat) to disk
    csvout = fullfile(pwd, conout, 'csv', [ edgew{edge} '.csv' ]);
    dlmwrite(csvout, omat(:,:,edgei(edge)), 'delimiter', ',');
    
    % copy label.json to each output - probably redundant
    copyfile('label.json', fullfile(pwd, conout));
    
    % write index.json w/ basic descriptors - what is this for?
    index = struct('filename', csvout, 'unit', edgew{edge}, 'name', edgew{edge}, 'desc', [ 'FiNE - ' edgew{edge} ' weighted connectivity' ]);
    savejson('', index, fullfile(pwd, conout, 'index.json'));
    % this is unqiqe per node modality (?) so 
    
    % make the new output directory for json graphs
    mkdir(fullfile(pwd, netout));

    % jsongraph filename to write to disk
    jsnout = fullfile(pwd, netout, 'network.json');
        
    % write jsongraph to disk
    fnCreateJsonGraph(netw, edgew{edge}, jsnout);
    gzip(jsnout); % gzip by defualt
    delete(jsnout); % remove the unzipped file

end

out = 'done.';

end
