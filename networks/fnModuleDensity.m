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
%                  'sum'  - adds up the values within a module - used with
%                           virtual lesion EMD output for module level statistic
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

if(~exist('flag', 'var') || isempty(flag))
    flag = 1;
end

% determine the number of modules
labels = unique(ci);
mnum = size(labels, 1);

% create empty output matrix
Mden = nan(mnum);

% grab every unique pair of modules
xy = nchoosek(1:mnum, 2);

% add the diagonal if it should be kept
if flag == 1
    pairs = sortrows([ [ 1:mnum; 1:mnum ]'; xy ]);
else
    pairs = xy;
end

clear xy

% only find the averages on observations in the upper diagonal
Mask = logical(triu(ones(size(M, 1)), 1));
mask = cell(size(pairs, 1), 1);

% for every module combination
for ii = 1:size(pairs, 1)
    
    % pull the module indices explicitly
    m1 = pairs(ii, 1);
    m2 = pairs(ii, 2);
    
    % pull subset of the variables for the module intersection
    iind = logical(ci == m1);
    jind = logical(ci == m2);
    
    % create the full matrix of upper / lower module values
    ind1 = iind * jind';
    ind2 = jind * iind';
    ind = ind1 + ind2;
    
    % combine upper diagonal with specific module
    mask{ii} = Mask & ind;
    
    % subset the data
    data = M(mask{ii});
    
    % define the type of edge summary to compute
    switch method
        
        case 'mean'
            Mden(m1, m2) = nanmean(data(:));
            Mden(m2, m1) = nanmean(data(:));
            
        case 'median'
            Mden(m1, m2) = median(data(:), 'omitnan');
            Mden(m2, m1) = median(data(:), 'omitnan');
            
        case 'std'
            Mden(m1, m2) = std(data(:), [], 'omitnan');
            Mden(m2, m1) = std(data(:), [], 'omitnan');
            
        case 'var'
            Mden(m1, m2) = var(data(:), [], 'omitnan');
            Mden(m2, m1) = var(data(:), [], 'omitnan');
            
        case 'boot'
            
            % define size and output of bootstrapping procedure
            nrep = 10000;
            nobs = size(data(:), 1);
            boot = zeros(nrep, 1);
            
            % for a fixed, reasonable number of permutations
            for perm = 1:nrep
                
                % sample data with replacement
                p = randsample(data, nobs, true);
                
                % store the mean of the random samples
                boot(perm) = nanmean(p);
                
            end
            
            % compute the mean of means for bootstrapped module value
            Mden(m1, m2) = nanmean(boot);
            Mden(m2, m1) = nanmean(boot);
            
        case 'sum'
            Mden(m1, m2) = nansum(data(:));
            Mden(m1, m2) = nansum(data(:));
            
        otherwise
            error('Invalid method of summarizing provided.');
    end
    
end

clear ii m1 m2 iind jind ind1 ind2 ind data rmean nobs boot

% error if a nan is found
if(any(isnan(Mden)))
    error('One or more averages is NaN. This is impossible. Inspect the input data.')
end

% compute a mask for the output
m1mask = eye(mnum); 
m2mask = ones(mnum) .* ~eye(mnum); 

% find where non-zero entries are of within / between modules
ff1 = find(m1mask);
ff2 = find(m2mask);

% mean property of within / between module blocks
% note: modules composed of a single node may have density 0
Mden_in = nanmean(Mden(ff1));
Mden_bw = nanmean(Mden(ff2));

% lowest within vs highest between
Mden_in_lo = min(nonzeros(Mden(ff1)));
Mden_bw_hi = max(Mden(ff2));

end
