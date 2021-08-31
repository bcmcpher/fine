function [ done ] = run_networks(subject, parc, life, clean)
%[ done ] = run_networks(subject);
% Just run what I need as a single job - BL is too camparmentalized to
%   usefully run this with data types.
%

%% fixed parameters

projdir = '/N/slate/bcmcpher/can-nets';
subjdir = fullfile(projdir, 'subjects');
subj = fullfile(subjdir, subject);

maxDist = 2;
minStrm = 5;
minLength = 10;
maxVars = 4;
maxLengthStd = 4;
numNodes = 100;
maxIter = 4;

%% load the data

disp('Loading brain data...')

% load the fe data
load(fullfile(subj, 'output_fe.mat'), 'fe');

% apparently vistasoft can't parse strings returned by fullfile as filenames (?)
% need to convert to strings to load files, but only in an interactive
%    session - this isn't necessary if it's called as a fxn (wtf...)
% convertStringsToChars(fullfile())

%track = dtiImportFibersMrtrix(fullfile(subj, 'tracks.tck'));
nodes = niftiRead(convertStringsToChars(fullfile(subj, 'parc.nii.gz')));

% load microstructure
fa = niftiRead(convertStringsToChars(fullfile(subj, 'fa.nii.gz')));
md = niftiRead(convertStringsToChars(fullfile(subj, 'md.nii.gz')));
odi = niftiRead(convertStringsToChars(fullfile(subj, 'odi.nii.gz')));
ndi = niftiRead(convertStringsToChars(fullfile(subj, 'ndi.nii.gz')));
isovf = niftiRead(convertStringsToChars(fullfile(subj, 'isovf.nii.gz')));

% load the passed labels
label = loadjson(fullfile(subj, 'label.json'));

% extract data from the fe structure

disp('Extracting Streamline data...');

% load fibers from fg
fg = feGet(fe, 'fibers acpc');

% pull the fiber weights
wght = feGet(fe, 'fiber weights');

% pull the model
M = feGet(fe, 'model');

% optionally subset to only positively weighted fibers
if life
    nonzero = wght > 0;
    nonzidx = find(nonzero); % because the sparse tensor can't be logically indexed
    fg.fibers = fg.fibers(nonzero);
    wght = wght(nonzero);
    M.Phi = M.Phi(:,:,nonzidx);
end

% pull the number of directions
nTheta = feGet(fe, 'nbvals');

% pull the demeaned diffusion signal
dsig = feGet(fe, 'dsigdemeaned by voxel');

% pull the mean diffusion signal
S0 = feGet(fe, 'b0signalimage');

% load the dictionary for angle computation
dict = feGet(fe, 'orient');

% load rmse volume if it exists, otherwise create it
if ~exist(fullfile(subj, 'rmse.nii.gz'), 'file')
    rmse = fnRmseVolume(fe);
    rmse.fname = 'rmse.nii.gz';
    niftiWrite(rmse, fullfile(subj, 'rmse.nii.gz'));
else
    rmse = niftiRead(convertStringsToChars(fullfile(subj, 'rmse.nii.gz')));
end

% clear fe for memory usage
clear fe nonzero nonzidx

% fix the label names

disp('Fixing labels...');

% extract clean the names and labels
names = cellfun(@(x) { x.name x.voxel_value }, label, 'UniformOutput', false);

% clean the names and labels
names = cellfun(@(x) {strrep(x{1}, '.label', ''), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_LH_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_RH_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '+', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '.', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_L_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_R_', '_'), x{2} }, names, 'UniformOutput', false);
names = cellfun(@(x) {strrep(x{1}, '_ROI', ''), x{2} }, names, 'UniformOutput', false);

% deal w/ maTT being wrong

% pull the labels - assume 0-max cortical
mlab = max(cellfun(@(x) x{2}, names));

% build the subcortical names
lhsub = {'lh_Thalamus'; 'lh_Caudate'; 'lh_Putamen'; 'lh_Pallidum'; 'lh_Hippocampus'; 'lh_Amygdala'; 'lh_Accumbens'};
rhsub = strrep(lhsub, 'lh_', 'rh_');
wbsub = [ lhsub; rhsub ];

% build the subcortical name/values
subnm = cellfun(@(x,y) {x, y}, wbsub, num2cell([mlab+1:mlab+14]'), 'UniformOutput', false);

% append the name/value combination
names = [ names'; subnm ];

clear lhsub rhsub subnm

% switch parc
%     case 'hcp'
%         load('/N/slate/bcmcpher/can-nets/data/hcp-names.mat', 'names');
%     otherwise
%         error('Invalid parcellation requested.');
% end

%% define the edges

% count
% - volume / fa / noddi / rmse
% - tract profile
% - vl
% - links
% cleaned
% - volume / fa / noddi / rmse
% - tract profile
% - vl
% - links

% create the network object
netw = fnCreateEdges(nodes, fg, names, maxDist, minStrm, 'weights', wght);

% optionally clean the edges
if clean
    netw = fnCleanEdges(netw, fg, minLength, maxVars, maxLengthStd, numNodes, maxIter, minStrm);
end

%% average edge properties

% estimate the central tendency of microstructure within each edge
netw = fnAveragePropertyEdges(netw, fg, {fa, md, odi, ndi, isovf, rmse}, {'fa', 'md', 'odi', 'ndi', 'isovf', 'rmse'}, true);

%% average profile properties

% estimate the average fa within each edge
netw = fnTractProfileEdges(netw, fg, fa, 'fa', numNodes, minStrm);

netw = fnTractProfileEdges(netw, fg, {fa, md, odi, ndi, isovf, rmse}, {'fa', 'md', 'odi', 'ndi', 'isovf', 'rmse'}, numNodes, minStrm);

% estimate the average noddi within each edge
netw = fnTractProfileEdges(netw, fg, nd, 'noddi', numNodes, false);

% estimate the average RMSE within each edge
netw = fnTractProfileEdges(netw, fg, rmse, 'rmse', numNodes, false);

%% profile shapes

netw = fnTractCurveEdges(netw, fg, numNodes, minNum);

%% link network

netw1 = fnCreateLinkNetwork(netw, 'dice', 'volume'); % volume
netw2 = fnCreateLinkNetwork(netw, 'mi', 'fa'); % MI of FA
netw3 = fnCreateLinkNetwork(netw, 'mi', 'nd'); % MI of NODDI
netw4 = fnCreateLinkNetwork(netw, 'angle', 'life', M.Phi, dict); % angle of intersection

%% virtual lesion

if life    
    % run virtual lesion if LiFE is incorporated
    netw = fnVirtualLesionEdges(netw, M, wght, dsig, nTheta, S0);
end

%% create outputs

% create matrix field for export of all computed edge weights
netw = fnComputeMatrixEdges(netw);

% create connectivity matrices
[ omat, olab ] = fnCreateAdjacencyMatrices(netw);

%% save outputs



end
