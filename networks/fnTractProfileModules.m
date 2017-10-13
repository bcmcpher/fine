function [ mnten, seten ] = fnTractProfileModules(tptens, srt, method)
%fnTractProfileModules() takes a tract profile tensor and averages the
% profiles into modules for a summary. Returns mean and standard error of
% the modules.
%   

% grab number of nodes from input
nnodes = size(tptens, 3);

% find maximum number of modules
nmod = max(srt);

% create empty output matrix
mnten = zeros(nmod, nmod, nnodes);
seten = zeros(nmod, nmod, nnodes);

% for every module 'row'
for i = 1:nmod
    
    % pull the 'row' indices for the module
    iind = logical(srt == i);
    
    % for every module 'col'
    for j = 1:nmod
        
        % pull the 'col' indices for the module
        jind = logical(srt == j);
        
        % subset the input matrix by modules
        Mij = tptens(iind, jind, :);
                
        % initialize output data
        moddat = nan(1, nnodes);
        iter = 1;
        
        % for every unique profile
        for k1 = 1:size(Mij, 1)
            for k2 = 1:size(Mij, 2)
                            
                % drop dimensions before adding it to data
                dat = squeeze(Mij(k1, k2, :))';
            
                % stack the data
                moddat = [ moddat; dat ];
                
                % iterate counter
                iter = iter + 1;
                
            end
        end
        
        % drop the initializing nan row
        moddat = moddat(2:end, :);
        
        % define the type of edge summary to compute
        switch method
            
            case 'mean'
                
                % directly compute mean and standard error
                mnten(i, j, :) = nanmean(moddat, 1);
                seten(i, j, :) = nanstd(moddat, 1) / size(moddat, 1);
                
            case 'boot'
                
                % define size and output of bootstrapping procedure
                nrep = 10000;
                nobs = size(moddat, 1);
                boot = zeros(nrep, nnodes);
                
                % for a fixed, reasonable number of permutations
                for perm = 1:nrep
                    
                    % sample data with replacement
                    p = randsample(1:nobs, nobs, true);
                    
                    % store the mean of the random samples
                    boot(perm, :) = nanmean(moddat(p, :));
                    
                end
                
                % compute the mean of means for bootstrapped module value
                mnten(i, j, :) = nanmean(boot, 1);
                seten(i, j, :) = nanstd(boot, 1) ./ nrep;
                
            otherwise
                
                error('Invalid method of summarizing provided.');
                
        end
    end
end

end


