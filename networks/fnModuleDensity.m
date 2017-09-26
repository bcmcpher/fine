function [ Mden, Mden_in, Mden_bw, Mden_in_lo, Mden_bw_hi ] = fnModuleDensity(M, ci, method, flag)
%fnModuleDensity() takes an input matrix and vector of node assignments and
% summarizes the individual edges within the assignments.
%
% INPUTS:
%     M      - input matrix to compute modules for
%     ci     - vector assigning nodes of M to modules
%     method - how edges are summarized within modules. 
%              Either:
%                  'mean' - the average of the edge weights in the modules
%                  'boot' - bootstrap the average edge weight in the module
%     flag   - whether or not to discount diagonal; assumes they're square
%
% OUTPUTS:
%     Mden       - the computed module density matrix
%     Mden_in    - mean density within module blocks
%     Mden_bw    - mean density between module blocks
%     Mden_in_lo - lowest within module mean
%     Mden_bw_hi - highest between module mean
%   
% TODO:
% - other 'method' options?
% - default arguments
%
% EXAMPLE:
%
% % load data
% parc = niftiRead('labels.nii.gz');
% fg        = feGet(fe, 'fibers acpc');
% fibers    = fg.fibers;
% fibLength = fefgGet(fg, 'length');
% weights   = feGet(fe, 'fiberweights');
% yeoROIs   = dlmread('yeoLabels.txt');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create network
% [ omat, olab ] = fnCreateAdjacencyMatrix(pconn, 'nzw');
%
% % compute module network
% [ Mden, Mden_in, Mden_bw, Mden_in_lo, Mden_bw_hi ] = fnModuleDensity(omat(:,:,1), yeoROIs, 'mean', 1);
%
% Olaf Sporns and Brent McPherson (c), 2017 - Indiana University
%

% determine the number of modules
mnum = max(ci);

% create empty output matrix
Mden = zeros(mnum);

% for every module 'row'
for i = 1:mnum
    
    % pull the 'row' indices for the module
    iind = logical(ci == i);
    
    % for every module 'col'
    for j = 1:mnum
        
        % pull the 'col' indices for the module
        jind = logical(ci == j);
        
        % subset the input matrix by modules
        Mij = M(iind, jind);
        
        % if its the same module and main diagonal is not discounted
        if ((i == j) && (flag == 0))
            
            % it's just the subset 
            Mij = M(iind, jind);
            
        end
        
        % if it's the same module and the main diagonal is discounted
        if ((i == j) && (flag == 1))
            
            % subset input matrix by modules
            Mij = M(iind, jind);
            
            % create and apply the mask to the diagonal
            mask = ones(size(Mij)) .* ~eye(size(Mij, 1));
            Mij = Mij(logical(mask));
            
        end
        
        % define the type of edge summary to compute
        switch method
            case 'mean'
                Mden(i, j) = nanmean(Mij(:));
            case 'boot'
                
                % linearize data
                data = Mij(:);
                
                % define size and output of bootstrapping procedure
                nrep = 10000;
                nobs = size(data, 1);
                boot = zeros(nrep, 1);
                
                % for a fixed, reasonable number of permutations
                for perm = 1:nrep
                    
                    % sample data with replacement
                    p = randsample(data, nobs, true);
                    
                    % store the mean of the random samples
                    boot(perm) = nanmean(p);
                    
                end
                
                % compute the mean of means for bootstrapped module value
                Mden(i, j) = nanmean(boot);
                
            otherwise
                error('Invalid method of summarizing provided.');
        end
    end
end

% replace missing with zeros
Mden(isnan(Mden)) = 0;

% compute a mask for the output
m1mask = eye(mnum); 
m2mask = ones(mnum) .* ~eye(mnum); 

% find where non-zero entries are of within / between modules
ff1 = find(m1mask);
ff2 = find(m2mask);

% mean density of within / between module blocks
% note: modules composed of a single node may have density 0
Mden_in = nanmean(Mden(ff1));
Mden_bw = nanmean(Mden(ff2));

% lowest within vs highest between
Mden_in_lo = min(nonzeros(Mden(ff1)));
Mden_bw_hi = max(Mden(ff2));

end
