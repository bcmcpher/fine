function [ vmat, olab, vlout, pcmodl ] = fnModuleVirtualLesion(fe, pconn, label, srt, norm)
%fnModuleVirtualLesion() perfurms virtual lesion on all connections within
% or between a network module(s).
%   

% parse optional arguments
if(~exist('norm', 'var') || isempty(norm))
    norm = 0;
    display('Computing raw virtual lesion.');
else
    norm = 1;
    display('Computing normalized virtual lesion.');
end

% compute the size of everything here
szsrt = size(srt, 1);
szpcn = size(pconn, 1);

% recreate pconn indices to grab labels correctly
pcindx = nchoosek(1:szsrt, 2);

% size of reconstructed pconn indices
szpci = size(pcindx, 1);

% error if the wrong number of labels is passed
if szpci ~= szpcn
    error('The number of nodes assigned to modules does not match the number of nodes in the network.');
end

% find maximum number of modules
nmod = max(srt);

% find lower diagonal (between modules)
lowd = nchoosek(1:nmod, 2);

% find diagonal (within modules)
diag = [ 1:nmod; 1:nmod ]';

% combine the module indices to compare
pairs = [ diag; lowd ];

% size of unique module combinations
szprs = size(pairs, 1);

% catch labels corresponding to module index
indx = repmat({struct('mod1', [], 'mod2', [], 'indices', [], 'lengths', [], 'weights', [])}, [szprs, 1]);

display(['Assigning all connections to one of ' num2str(szprs) ' unique modules...']);

% preallocate connection assignment array 
pcmodl = zeros(szpci, 1);

% grab pairwise module assignment of every connection
parfor conn = 1:szpci
    
    % grab the module labels
    roi1 = srt(pcindx(conn, 1));
    roi2 = srt(pcindx(conn, 2));
    
    % grab and return the module this connection is in
    % this will try both orders of ROI indexing (only 1 will work)
    try
        [ ~, pcmodl(conn) ] = intersect(pairs, [ roi1 roi2 ], 'rows');
    catch
        [ ~, pcmodl(conn) ] = intersect(pairs, [ roi2 roi1 ], 'rows');
    end
    
end

display(['Condensing streamline indices from ' num2str(size(pconn, 1)) ' connections into ' num2str(size(indx, 1)) ' modules...']);

% for every module
parfor mod = 1:szprs
    
    % for every connection in every module
    modi = find(pcmodl == mod);
    
    % grab all the streamline indices
    for conn = 1:size(modi)
        
        % grab the indices
        indx{mod}.mod1 = pairs(mod, 1);
        indx{mod}.mod2 = pairs(mod, 2);
        indx{mod}.indices = [ indx{mod}.indices; pconn{modi(conn)}.(label).indices ];
        indx{mod}.lengths = [ indx{mod}.lengths; pconn{modi(conn)}.(label).lengths ];
        indx{mod}.weights = [ indx{mod}.weights; pconn{modi(conn)}.(label).weights ];
        
    end
    
end

display(['Computing virtual lesion for every one of ' num2str(size(indx, 1)) ' network modules...']);

% preallocate output
vmat = zeros(nmod, nmod, 4);
vlout = cell(szprs, 1);

parfor vl = 1:szprs
    
    % compute virtual lesion as either raw or normalized
    switch norm
        case 0
            [ wVL, woVL ] = feComputeVirtualLesion(fe, indx{vl}.indices);
            vlout{vl} = feComputeEvidence(woVL, wVL);
        case 1
            [ wVL, woVL ] = feComputeVirtualLesion_norm(fe, indx{vl}.indices);
            vlout{vl} = feComputeEvidence_norm(woVL, wVL);
        otherwise
            error('Bad VL normalization call that can''t happen?');
    end
    
end

% collect virtual lesions into output matrix
for vl = 1:szprs
    
    % create upper diagonal output matrix
    vmat(pairs(vl, 1), pairs(vl, 2), 1) = vlout{vl}.em.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 2) = vlout{vl}.s.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 3) = vlout{vl}.kl.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 4) = vlout{vl}.j.mean;

    % create lower diagonal output matrix
    vmat(pairs(vl, 2), pairs(vl, 1), 1) = vlout{vl}.em.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 2) = vlout{vl}.s.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 3) = vlout{vl}.kl.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 4) = vlout{vl}.j.mean;

end

% output labels that do not change
olab = {'EMD', 'SOE', 'KDL', 'JfD'};

end

