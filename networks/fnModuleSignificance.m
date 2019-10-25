function [ Mden, fh, Mprm ] = fnModuleSignificance(M, ci, grp, Nperm, flag)
%fnModuleDensity() takes an input matrix and vector of node assignments and
% summarizes the individual edges within the assignments.
%
% INPUTS:
%     M      - input matrix to compute modules for
%     ci     - vector assigning nodes of M to modules
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
% assumes the modules are labeled 1-N with no gaps
labels = unique(ci);
mnum = size(labels, 1);

% LABELS NEED TO BE ASSIGNED A GROUP - BRAIN OR BEHAVIOR
% RESAMPLING FOR NULL MODULE DISTRIBUTIONS NEEDS TO BE RESAMPLED WITHIN
% GROUP: BRAIN-BRAIN, BEHAVIOR-BEHAVIOR, BRAIN-BEHAVIOR

% create empty output matrix
Mden = nan(mnum);
Mprm = nan(mnum, mnum, Nperm);

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
nsize = nan(size(pairs, 1), 1);
grps = nan(size(pairs, 1), 1);
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
    
    % grab the null size
    nsize(ii) = size(data, 1);
    
    % pull the true mean
    rmean = nanmean(data);

    % define size and output of bootstrapping procedure
    nobs = size(data, 1);
    boot = zeros(Nperm, 1);
    
    % for a fixed, reasonable number of permutations
    for perm = 1:Nperm
        
        % sample data with replacement
        p = randsample(data, nobs, true);
        
        % store the mean of the random samples
        boot(perm) = nanmean(p);
        
    end
    
    % assign the real value out
    Mden(m1, m2) = rmean;
    Mden(m2, m1) = rmean;
    
    % assign the bootstrap out
    Mprm(pairs(ii, 1), pairs(ii, 2), :) = boot;
    Mprm(pairs(ii, 2), pairs(ii, 1), :) = boot;
    
    % assign sum as unique groups
    grps(ii) = grp(m1) + grp(m2);
    
end

grps = grps - 1;

clear ii m1 m2 iind jind ind1 ind2 ind data rmean nobs boot

% error if a nan is found
if(any(isnan(Mden)))
    error('One or more averages is NaN. This is impossible. Inspect the input data.')
end

%% plot the null distributions with the real values

nvals = size(M, 1);
%ngrps = size(unique(grps), 1);

%wnmsk = cell(size(unique(grps), 1), 1);
wnmsk = zeros(nvals, nvals, size(unique(grps), 1));

for ii = 1:size(pairs, 1)
    wnmsk(:, :, grps(ii)) = wnmsk(:, :, grps(ii)) + mask{ii};
end

% alternative / additional?
% get the null distribution of the global mean by pulling random data into modules of the same real size
% plot the true mean of each module to determine where they individually
% fall on the tails
ngm = nan(Nperm, size(nsize, 1));
for ii = 1:Nperm
    for jj = 1:size(nsize, 1)
        tdat = randsample(M(mask{jj}), nsize(jj), 'true');
        ngm(ii, jj) = nanmean(tdat);        
    end
end

% build the indices into the subplots to match the diagonal
pidx = find(tril(ones(mnum)));
%mnmn = mean(Mden(:));

fh = figure;

% for every module
for ii = 1:size(pairs, 1)
    
    m1 = pairs(ii, 1);
    m2 = pairs(ii, 2);
    
    % subplot the histogram of the null averages w/ the real value in red
    % null is estimated w/in domain, controlling for module size
    subplot(mnum, mnum, pidx(ii));
    hist(squeeze(ngm(:, 1)), 32);
    line([ Mden(m1, m2) Mden(m1, m2) ], [ 0 100 ], 'color', 'red');
    set(gca, 'XTick', [], 'YTick', []);
    
end

figure;
for ii = 1:size(unique(grps), 1)
    
    subplot(1, size(unique(grps), 1), ii);
    
    % pull the right split of data
    idx = grps == ii;
    xy = pairs(idx, :);
    
    % subset the null data
    tdat = ngm(:, idx);
    
    % histogram of the within domain null distribution
    hist(tdat(:), 64);
    
    % for every observation within the domain
    for jj = 1:size(xy, 1)
        % plot the line of the true mean
        line([ Mden(xy(jj, 1), xy(jj, 2)) Mden(xy(jj, 1), xy(jj, 2)) ], [ 0 3000 ], 'color', 'red');
    end
    
end

end
