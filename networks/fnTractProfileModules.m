function [ mnten, seten, mnmat, semat ] = fnTractProfileModules(tptens, srt, method, nrep)
%fnTractProfileModules() takes a tract profile tensor and averages the
% profiles into modules for a summary. Returns mean and standard error of
% the modules.
% 

% parse optional arguments
if(~exist('method', 'var') || isempty(method))
    method = 'mean';
end

% assign bootstrap replications by default
if(~exist('nrep', 'var') || isempty(nrep))
    nrep = 10000;
end

% grab number of nodes from input
nnodes = size(tptens, 3);

% find maximum number of modules
nmod = length(unique(srt));

% create empty output matrix
mnten = nan(nmod, nmod, nnodes);
seten = nan(nmod, nmod, nnodes);

% find lower diagonal (between modules)
uppd = nchoosek(1:nmod, 2);

% find diagonal (within modules)
diag = [ 1:nmod; 1:nmod ]';

% combine the module indices to compare
pairs = [ diag; uppd ];

% sort into upper diagonal indices
pairs = sortrows(pairs, [ 1 2 ]);

% create pairs/nmod x nnodes average matrices
mnmat = nan(size(pairs, 1), nnodes);
semat = nan(size(pairs, 1), nnodes);

% for every combination of modules, including diagonal
for ii = 1:size(pairs, 1)
    
    % pull module references
    mod1 = pairs(ii, 1);
    mod2 = pairs(ii, 2);
    
    % grab the logical indices for the module combination
    iind = logical(srt == mod1);
    jind = logical(srt == mod2);
    
    % subset the input matrix by modules
    Mij = tptens(iind, jind, :);

    % pull the indices of the matrix (not necessarily symmetric)
    [ x1, y1 ] = ind2sub(size(Mij, [ 1 2 ]), 1:prod(size(Mij,[ 1 2 ])));
    tpPairs = [ x1', y1' ];
    
    % reshape into profile in module x 100 matrix
    Mpf = nan(size(tpPairs, 1), nnodes);
    for jj = 1:size(tpPairs, 1)
        Mpf(jj, :) = squeeze(Mij(tpPairs(jj, 1), tpPairs(jj, 2), :))';
    end
    
    % drop all missing and all zeros
    Mpf(isnan(Mpf)) = 0; % set NaN to zero
    Mpf = Mpf(any(Mpf, 2), :); % drop all zeros
    
    switch method
        
        case {'mean', 'mn-sd'}
            
            % compute the mean / sd of the profiles
            tpmn = mean(Mpf, 'omitnan');
            tpse = std(Mpf, 'omitnan');
            
            % add to output
            mnten(mod1, mod2, :) = tpmn;
            seten(mod1, mod2, :) = tpse;
            mnmat(ii, :) = tpmn;
            semat(ii, :) = tpse;
            
            % if it's not on the diagonal
            if mod1 ~= mod2
                
                % add flipped 
                mnten(mod2, mod1, :) = fliplr(tpmn);
                seten(mod2, mod1, :) = fliplr(tpse);
                
            end
            
        case 'mn-se'
            
            % compute the mean / se of the profiles
            tpmn = mean(Mpf, 'omitnan');
            tpse = std(Mpf, 'omitnan') ./ sqrt(size(Mpf, 1));
            
            % add to output
            mnten(mod1, mod2, :) = tpmn;
            seten(mod1, mod2, :) = tpse;
            mnmat(ii, :) = tpmn;
            semat(ii, :) = tpse;
            
            % if it's not on the diagonal
            if mod1 ~= mod2
                
                % add flipped 
                mnten(mod2, mod1, :) = fliplr(tpmn);
                seten(mod2, mod1, :) = fliplr(tpse);
                
            end            
            
        case 'boot'
            
            % define size and output of bootstrapping procedure
            nobs = size(Mpf, 1);
            boot = zeros(nrep, nnodes);
            
            % for a fixed, reasonable number of permutations
            for perm = 1:nrep
                
                % sample data with replacement
                p = randsample(1:nobs, nobs, true);
                
                % store the mean of the random samples
                boot(perm, :) = mean(Mpf(p, :), 1);
                
            end
            
            % compute the mean of means for bootstrapped module value
            tpmn = mean(boot, 'omitnan');
            tpse = std(boot, 'omitnan') ./ nrep;
            
            % add to output
            mnten(mod1, mod2, :) = tpmn;
            seten(mod1, mod2, :) = tpse;
            mnmat(ii, :) = tpmn;
            semat(ii, :) = tpse;
            
            % if it's not on the diagonal
            if mod1 ~= mod2
                
                % add flipped 
                mnten(mod2, mod1, :) = fliplr(tpmn);
                seten(mod2, mod1, :) = fliplr(tpse);
                
            end            
                
        otherwise
            
            error('Invalid method of summarizing provided.');
            
    end
    
end
   
% set any final nan to zero
mnten(isnan(mnten)) = 0;
seten(isnan(seten)) = 0;
mnmat(isnan(mnmat)) = 0;
semat(isnan(semat)) = 0;
