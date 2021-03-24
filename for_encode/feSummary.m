function [ out ] = feSummary(fe, class)
%out = feSummary() Provides the summary output Dan wrote for an fe structure.
% Will optionally parse and return classified fiber breakdown as well. 
%
%   Return some summary values of an fe structure, including:
%       - normalized RMSE, RMSE
%       - total streamline length
%       - proportion of weighted / assigned streamlines
% 
%   Optionally, if a classification structure is passed the same relevant
%   numbers are generated for the classified streamlines.
%
% INPUTS:
%     fe    - an evaluated fiber evaluation structure (LiFE has been performed)
%     class - a fiber classification structure
% 
% OUTPUTS:
%     out   - a struct containing the summary statistics on the fe and
%             optional classification structure.
%
% Dan Bullock and Brent McPherson (c) 2018 Indiana University
%

if(~exist('class', 'var') || isempty(class))
    class = [];
    out = [];
end

disp('Extracting data...');

% get acpc streamlines
wb_fg = feGet(fe, 'fibers acpc');

% pull weights
weights = feGet(fe, 'fiberweights');

% grab the number of streamlines
nsl = size(weights, 1);

% gets positively weighted streamlines and their indexes
nzw = find(weights > 0);

% a test for the nan
nanNum = sum(isnan(weights));
if nanNum > 0
    fprintf('\n %i NaN weights detected.', nanNum)
end

% drop nan weighted fibers
nz_fg = wb_fg;
nz_fg.fibers = wb_fg.fibers(nzw);

%% compute fe structure summary

disp('Computing streamline lengths...');

% preallocate length output
wb_fg_lengths = nan(nsl, 1);

% computes length of all streamlines
for sl = 1:nsl
    wb_fg_lengths(sl) = sum(sqrt(sum(diff(wb_fg.fibers{sl}, 1, 2).^2)));
end

% length of positively weighted streamlines
nz_fg_lengths = wb_fg_lengths(nzw);

% grab fe name
out.name = fe.name;

disp('Computing global RMSE signal explained...');

% grab normalized RMSE
nrm_rmse = feGet(fe, 'voxrmses0norm');

% whole brain normed RMSE
out.rmse.nrm_total = sum(nrm_rmse, 'omitnan');

% grab raw RMSE
out.rmse.raw_total = feGet(fe, 'rmsetotal');

% additional descriptives of normed RMSE
out.rmse.nrm_vx_mn = mean(nrm_rmse, 'omitnan');
out.rmse.nrm_vx_sd = std(nrm_rmse, 'omitnan');
%out.rmse.nrm.vx_mx = max(nrm_rmse);
%out.rmse.nrm.vx_mn = min(nrm_rmse);
%out.rmse.nrm.vx_ct = length(nrm_rmse);

disp('Computing streamline descriptives...');

% compute total length of pre- and post-LiFE streamline connectome
out.length.wb_total = sum(wb_fg_lengths);
out.length.nz_total = sum(nz_fg_lengths);
out.length.proportion = out.length.wb_total / out.length.nz_total;

% compute proportion of nzw streamlines
out.life.wb_count = nsl;
out.life.nz_count = length(nzw);
out.life.proportion = out.life.nz_count / out.life.wb_count;

% mean, std, and max of streamline lengths
out.streamline.wb_mn_length = mean(wb_fg_lengths, 'omitnan');
out.streamline.wb_sd_length = std(wb_fg_lengths, 'omitnan');
out.streamline.wb_mx_length = max(wb_fg_lengths);

out.streamline.nz_mn_length = mean(nz_fg_lengths, 'omitnan');
out.streamline.nz_sd_length = std(nz_fg_lengths, 'omitnan');
out.streamline.nz_mx_length = max(nz_fg_lengths);

if ~isempty(class)
    
    % grab all classified indices
    indx = find(class.index);
    
    if length(class.index) == out.life.nz_count
        
        % if the number classified matches the nz number of streamlines
        fprintf ('\n Number of items in the classification structure suggests that classification was performed on pruned connectome.')
        fprintf ('\n Proceeding with this assumption...')
        
        % an unknown number of streamlines was input
        out.classified.wb_count = [];
        out.classified.wb_proportion = [];
        out.classified.wb_mn_length = [];
        out.classified.wb_sd_length = [];
        
        % grab the nz classified count
        out.classified.nz_count = length(find(class.index));
        
        % grab the nz classified proportion
        out.classified.nz_proportion = out.classified.nz_count / out.life.nz_count;

        % grab the nz classified length descriptives
        out.classified.nz_mn_length = mean(out.length.nz_total(indx), 'omitnan');
        out.classified.nz_sd_length = std(out.length.nz_total(indx), 'omitnan');
    
    % if the number classified matches the wb number of streamlines    
    elseif length(class.index) == out.life.wb_count        
        
        % grab the unclassified count and proportion
        out.classified.wb_count = length(indx);
        out.classified.wb_proportion = out.classified.wb_count / out.life.wb_count;
        
        % grab the unclassified streamline length descriptives
        out.classified.wb_mn_length = mean(wb_fg_lengths(indx));
        out.classified.wb_sd_length = std(wb_fg_lengths(indx));
        
        % find the intersection of the classified indices and weighted streamlines
        srv = indx(ismember(indx, nzw));
        
        % grab the nz classified count 
        out.classified.nz_count = length(srv);
        
        % grab the nz classified proportion
        out.classified.nz_proportion = out.classified.nz_count / out.life.nz_count;
        
        % grab the nz classified length descriptives
        out.classified.nz_mn_length = mean(nz_fg_lengths(srv), 'omitnan');
        out.classified.nz_sd_length = std(nz_fg_lengths(srv), 'omitnan');
        
    else
        
        fprintf('\n Number of streamlines classified is neither equal to total number of streamlines in wb fg nor total number of validated streamlines.')
        fprintf('\n Length of ''class.index'' should be exactly equal to the number of streamlines from input tractography structure.');
        
    end
    
end

end

