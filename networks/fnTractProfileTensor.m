function [ pmat, tpmat ] = fnTractProfileTensor(netw, prof, srt)
%fnCreateTractProfileTensor() creates an NxNxNodes tensor of all
% precomputed tract profiles. Empty profiles are composed of NaNs. 
%   
% INPUTS:
%     netw  - is the paired connections object to create adjacency matrices from.
%
%     prof  - a string indicating the profile to extract and combine into
%             the tensor. This is either the 'mslab' argument passed to
%             feTracProfilePairedConnections() or if a dt6 was used,
%             'dt6.<value>' where value is any of the following: fa, md, rd, ad, cl, cp, cs
%     
%     srt   - sort the indices 
%
% OUTPUTS:
%     pmat  - 3d array containing (nodes x nodes x profile) of tract profiles
%     tpmat - 2d array containing (connection x profile) tract profiles ordered by the upper diagonal 
% 
%
% EXAMPLE:
%
% % load data
% parc          = niftiRead('labels.nii.gz');
% fg            = feGet(fe, 'fibers acpc');
% fibers        = fg.fibers;
% fibLength     = fefgGet(fg, 'length');
% weights       = feGet(fe, 'fiberweights');
% dt            = dtiloadDt6('dt6.mat');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create a profile for every connection with non-zero weighted fibers
% pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', favol, 'fa');
%
% % create a profile for every connection with non-zero weighted fibers
% pconn = feTractProfilePairedConnections(fg, pconn, 'nzw', dt, 'dt6');
%
% % create tract profile tensor of FA
% pmat = feTractProfileTensor(pconn, 'nzw', 'fa');
%
% create tract profile tensor of dt6 FA
% pmat = feTractProfileTensor(pconn, 'nzw', 'dt6.fa');
%
% Brent McPherson (c), 2017 - Indiana University
%

%% parse arguments

if(~exist('prof', 'var') || isempty(prof))
    error('You have to request a kind of profile.');
end

if ~isfield(netw.edges{1}.profile, prof)
    error('The requested profile does not exist.');
end

if(~exist('srt', 'var') || isempty(srt))
    srt = [];
end

%% extrac the data from the network object

% look at first entry for the number of nodes
nnodes = size(netw.edges{1}.profile.(prof), 1);

% pull the number of nodes and pairs to index through
nlabs = size(netw.nodes, 1);
pairs = netw.parc.pairs;

% pre-allocate output to fill in profiles
pmat = zeros(nlabs, nlabs, nnodes);
tpmat = nan(size(pairs, 1), nnodes);

%% pull every profile

disp(['Extracting tract profiles from ''' prof '''...']);

% for every possible connection
for ii = 1:size(pairs, 1)
    
    % simple indices of unique edges
    grp1 = pairs(ii, 1);
    grp2 = pairs(ii, 2);
   
    % pull each profile
    tmp = netw.edges{ii}.profile.(prof);
    
    % if there is no tract fill, in profile as missings
    if isempty(tmp)
        tmp = nan(nnodes, 1);
    end
    
    % assign profile to profile tensor
    pmat(grp1, grp2, :) = tmp;
    pmat(grp2, grp1, :) = flipud(tmp);
    tpmat(ii, :) = tmp;
    
end

% sort the cortical nodes given a vector
if ~isempty(srt)
    pmat = pmat(srt, srt, :);
end

end
