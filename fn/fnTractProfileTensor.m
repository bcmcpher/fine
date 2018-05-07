function [ pmat, tpmat ] = fnTractProfileTensor(pconn, label, prof, srt)
%fnCreateTractProfileTensor() creates an NxNxNodes tensor of all
% precomputed tract profiles. Empty profiles are composed of NaNs. 
%   
% INPUTS:
%     pconn - is the paired connections object to create adjacency matrices from.
%
%     label - string indicating the fiber groups for which to create virtual lesions
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%             Additionally, this can be run after cleaning, resulting in
%             valid calls of 'all_clean' and 'nzw_clean', respectively.
%
%     prof  - a string indicating the profile to extract and combine into
%             the tensor. This is either the 'mslab' argument passed to
%             feTracProfilePairedConnections() or if a dt6 was used,
%             'dt6.<value>' where value is any of the following: fa, md, rd, ad, cl, cp, cs
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

% parse optional arguments
%if(~exist('nnodes', 'var') || isempty(nnodes))
%    nnodes = 100;
%end

if(~exist('srt', 'var') || isempty(srt))
    srt = [];
end

% look at first entry for the number of nodes
nnodes = size(eval([ 'pconn{1}.(label).profile.' prof ]), 1);

% find the number of unique labels
uniquelabels = zeros(size(pconn, 1), 1);
for ii = 1:size(pconn, 1)
    uniquelabels(ii, 1) = pconn{ii}.roi1;
end
nlabs = size(unique(uniquelabels), 1) + 1;

clear uniquelabels ii 

% re-create pairs based on number of nodes
pairs = nchoosek(1:nlabs, 2);

% pre-allocate output to fill in profiles
pmat = zeros(nlabs, nlabs, nnodes);
tpmat = nan(size(pairs, 1), nnodes);

display(['Extracting tract profiles from ' label ' ...']);

% for every possible connection
for ii = 1:size(pairs, 1)
    
    % simple indices of unique edges
    grp1 = pairs(ii, 1);
    grp2 = pairs(ii, 2);
   
    % pull each profile
    tmp = eval([ 'pconn{ii}.(label).profile.' prof ]);
    
    % if there is no tract fill, in profile as missings
    if isempty(tmp)
        tmp = nan(nnodes, 1);
    end
    
    % assign profile to profile tensor
    pmat(grp1, grp2, :) = tmp;
    pmat(grp2, grp1, :) = fliplr(tmp);
    tpmat(ii, :) = tmp;
    
end

% sort the cortical nodes given a vector
if ~isempty(srt)
    pmat = pmat(srt, srt, :);
end

end

