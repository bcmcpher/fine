%% create link network
% Brent McPherson
% 20170212
% 

%% load data

%load('cleaned_networks_20170130.mat', 'pconn');

%% generate link network

% link nodes
lnodes = length(pconn);

% create indices for every pair of connections
pairs = nchoosek(1:lnodes, 2);

% build an update list based on percentiles through length of pairs

% build empty cells to parallelize computation
ecell = zeros(lnodes, 1);
ccell = zeros(lnodes, 1);

% build empty matrices
emat = zeros(lnodes, lnodes);
cmat = zeros(lnodes, lnodes);

%parpool;

% for every connections overlap
parfor ii = 1:length(pairs)
    
    % simple index for each 
    ti1 = pairs(ii, 1);
    ti2 = pairs(ii, 2);
    
    % find size of intersection - denominator of Dice coeff
    esize = size(unique([pconn{ti1}.vxend; pconn{ti2}.vxend], 'rows'), 1);
    csize = size(unique([pconn{ti1}.vxcln; pconn{ti2}.vxcln], 'rows'), 1);
    %esize = size(pconn{ti1}.vxend, 1) + size(pconn{ti2}.vxend, 1);
    %csize = size(pconn{ti1}.vxcln, 1) + size(pconn{ti2}.vxcln, 1);
            
    % find the intersection - numerator of Dice coeff
    eint = ismember(pconn{ti1}.vxend, pconn{ti2}.vxend, 'rows');
    cint = ismember(pconn{ti1}.vxcln, pconn{ti2}.vxcln, 'rows');
    
    % make numerator a single number
    enum = sum(eint);
    cnum = sum(cint);
    
    % build Dice coeff values for assignment
    ecell(ii) = (2 * enum) / esize;
    ccell(ii) = (2 * cnum) / csize;

end

clear ii ti1 ti2 esize csize eint cint enum cnum

%% assemble matrix

for ii = 1:length(pairs)
    
    % simple index for each 
    ti1 = pairs(ii, 1);
    ti2 = pairs(ii, 2);
    
    % build matrices
    emat(ti1, ti2) = ecell(ii);
    emat(ti2, ti1) = ecell(ii);
    
    cmat(ti1, ti2) = ccell(ii);
    cmat(ti2, ti1) = ccell(ii);

end

clear ii ti1 ti2 

% fix nan/inf values to zero
emat(isinf(emat)) = 0;
emat(isnan(emat)) = 0;

cmat(isinf(cmat)) = 0;
cmat(isnan(cmat)) = 0;

%% compute stats on link matrix

% compute reduced stats due to increased matrix size
[ efns, efdat ] = feLinkNetworkStats(emat, 50, 1);
[ cfns, cfdat ] = feLinkNetworkStats(cmat, 50, 1);

%% simple plot    

plt = log10(emat);
plt = cfdat.agree;

fh = figure();
colormap('hot');
imagesc(plt);
axis('square'); axis('equal'); axis('tight');
colorbar;
set(gca, 'XTickLabel', '', 'YTickLabel', '', 'XTick', [], 'YTick', []);

%caxis(cmap);
%line([34.5 34.5], [0.5 68.5], 'Color', [0 0 1]);
%line([0.5 68.5], [34.5 34.5], 'Color', [0 0 1]);
%line([68.5 0.5], [68.5 0.5], 'Color', [0 0 1]);
    
%% cosine similarity metric development

% % the data values
% x = pconn{2}.vxend;
% y = pconn{15}.vxend;
% 
% % manually?
% theta = (x * y') / (norm(x, 2) * norm(y, 2));
% theta = acos((x * y') / (norm(x, 2) * norm(y, 2)));
% 
% % with built-in?
% theta = pdist([x; y], 'cosine'); % what do I do with the vector?
