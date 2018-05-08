function [ vmat, olab, indx, pcmodl ] = fnModuleVirtualLesion(fe, pconn, label, srt, norm)
%fnModuleVirtualLesion() perfurms virtual lesion on all connections within
% or between a network module(s).
%   
% INPUTS:
%     fe    - a fit LiFE model
%     pconn - a paired connection structure assigning every streamline 
%             to a connection
%     label - which set of indices to assign from pconn. Either
%               - all: all streamlines assigned
%               - nzw: non-zero weighted
%               - zwr: zero weighted / removed (cannot virtual lesion)
%     srt   - node community assignment, size(pconn, 1) x 1 array
%             determines the module community of every connection
%     norm  - (optional, default = 0) 1 or 0, whether or not to normalize 
%             the virtual lesion computation
%
% OUTPUTS:
%    vmat   - a nmod x nmod x 4 array of the requested modules with each
%             virtual lesion measure 
%    olab   - a 1 x 4 string array indicating the measure along the 4th
%             dimension of vmat
%    indx   - the raw data used within the function to compute the vl
%    pcmodl - an array with the module assignment of every connection in pconn
%
% Brent McPherson, (c) 2017, Indiana University
%

% parse optional vl normalization arguments
if(~exist('norm', 'var') || isempty(norm))
    norm = 0;
    disp('Computing raw virtual lesion.');
else
    norm = 1;
    disp('Computing normalized virtual lesion.');
end

% compute the size of node assignments and the number of connections
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
nmod = size(unique(srt), 1);

% find upper diagonal (between modules)
uprd = nchoosek(1:nmod, 2);

% find diagonal (within modules)
diag = [ 1:nmod; 1:nmod ]';

% combine the module indices to compare
pairs = [ diag; uprd ];

% size of unique module combinations
szprs = size(pairs, 1);

disp(['Assigning all connections to one of ' num2str(szprs) ' unique modules...']);

% preallocate connection assignment array 
pcmodl = nan(szpci, 1);

% catch labels corresponding to module index
indx = repmat({struct('mod1', [], 'mod2', [], 'indices', [], 'lengths', [], 'weights', [])}, [szprs, 1]);

% grab module assignment of every connection
for conn = 1:szpci
    
    % grab the module labels - srt value
    roi1 = srt(pcindx(conn, 1));
    roi2 = srt(pcindx(conn, 2));
    
    % grab and return the module this connection is in
    % this will try both orders of ROI indexing (only 1 will work)
    try
        % use the row index in pairs to assign the community label
        [ ~, pcmodl(conn) ] = intersect(pairs, [ roi1 roi2 ], 'rows');
    catch
        [ ~, pcmodl(conn) ] = intersect(pairs, [ roi2 roi1 ], 'rows');
    end
    
    % grab module index
    mod = pcmodl(conn);
    
    % assign module pairs - is reassigned every loop...
    indx{mod}.mod1 = pairs(mod, 1);
    indx{mod}.mod2 = pairs(mod, 2);
    
    % grab the indices / lengths / weights for every streamline in each connection
    indx{mod}.indices = [ indx{mod}.indices; pconn{conn}.(label).indices ];
    indx{mod}.lengths = [ indx{mod}.lengths; pconn{conn}.(label).lengths ];
    indx{mod}.weights = [ indx{mod}.weights; pconn{conn}.(label).weights ];
    
end

clear conn roi1 roi2 mod

% total summaries
tcon = 0;
tfib = zeros(szprs, 1);

% report the number of connections and total streamlines of each module
for mod = 1:szprs
    
    % pull values from indx
    nconn = size(find(pcmodl == mod), 1);
    nfib = size(indx{mod}.indices, 1);
    
    % print out each module
    disp([ 'Module ' num2str(mod) ' has ' num2str(nconn) ' connections and ' num2str(nfib) ' streamlines.' ]);
    
    % iterate total counts
    tcon = tcon + nconn;
    tfib(mod) = nfib; 
    
end

% report total connections and streamlines assigned (should be all of them)
disp([ 'A total of ' num2str(tcon) ' connections and ' num2str(sum(tfib)) ' streamlines were assigned.' ]);  

clear tcon tfib mod nconn nfib

disp(['Computing virtual lesion for every one of ' num2str(size(indx, 1)) ' network modules...']);

parfor vl = 1:szprs
    
    if all(indx{vl}.weights == 0)
        
        warning('Connection %d has no positive weighted streamlines. VL not computed.', vl)
        indx{vl}.vl.em.mean = 0;
        indx{vl}.vl.s.mean = 0;
        indx{vl}.vl.kl.mean = 0;
        indx{vl}.vl.j.mean = 0;
        
        continue
    
    end
        
    % compute virtual lesion as either raw or normalized
    switch norm
        case 0
            [ wVL, woVL ] = feComputeVirtualLesion(fe, indx{vl}.indices);
            indx{vl}.vl = feComputeEvidence(woVL, wVL);
        case 1
            [ wVL, woVL ] = feComputeVirtualLesion_norm(fe, indx{vl}.indices);
            indx{vl}.vl = feComputeEvidence_norm(woVL, wVL);
        otherwise
            error('Bad VL normalization call that can''t happen?');
    end
    
end

% preallocate output
vmat = zeros(nmod, nmod, 4);

% collect virtual lesions into output matrix
for vl = 1:szprs
    
    % create upper diagonal output matrix
    vmat(pairs(vl, 1), pairs(vl, 2), 1) = indx{vl}.vl.em.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 2) = indx{vl}.vl.s.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 3) = indx{vl}.vl.kl.mean;
    vmat(pairs(vl, 1), pairs(vl, 2), 4) = indx{vl}.vl.j.mean;

    % create lower diagonal output matrix
    vmat(pairs(vl, 2), pairs(vl, 1), 1) = indx{vl}.vl.em.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 2) = indx{vl}.vl.s.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 3) = indx{vl}.vl.kl.mean;
    vmat(pairs(vl, 2), pairs(vl, 1), 4) = indx{vl}.vl.j.mean;

end

% output labels that do not change
olab = {'EMD', 'SOE', 'KLD', 'JfD'};

end

