function [ netw, glob, node, nets ] = fnNetworkStatistics(netw, nperm, iter)
%[ netw, glob, node, nets ] = fnNetworkStatistics(netw, nperm, iter);
%   This estimates a series of global, nodal, and edge-wise network
%   statistics. It will optionally randomly sort the network for estimating
%   the significance of some properties. 
%   
%   Some properties are preserved through the random sorting, so a pvalue
%   above threshold for some properties is normal and not indicative an
%   issue. 
%   
%   TODO:
%   - select how the random sorting happens? maybe varies based on fxn
%   - external fxns for calling different combinations of properties?
%   - handle previously stored values when they're encountered (clobber)
%   - select an edge weight for doing the weighted estimates on a specific value
%

%% parse default arguments

if(~exist('nperm', 'var') || isempty(nperm))
    nperm = [];
end

if(~exist('iter', 'var') || isempty(iter))
    iter = 10;
end

%% create data

% pull some basic shape info
nnodes = size(netw.nodes, 1);
nedges = size(netw.edges, 1);

% create fast, simple matrix
mat = squareform(cellfun(@(x) length(x.fibers.indices), netw.edges));

% binarize it
bin = weight_conversion(mat, 'binarize');

%% run null permutations if requested

if ~isempty(nperm)
    
    disp('Building null matrices...');
    
    % preallocate null array
    null = cell(nperm, 1);
    
    for perm = 1:nperm
        
        % randomly permute matrix
        rmat = randmio_und_connected(bin, iter);
        % it is unclear what properties would be preserved in this shuffle
        
        % store permutations
        null{perm} = rmat;
        
    end 
end

% run stats on network
disp('Estimating the observed network statistics...');
[ glob, node, nets ] = runStatistics(bin);

% pull the names of measures in statistics objects
gmeas = fieldnames(glob);
nmeas = fieldnames(node);
emeas = fieldnames(nets);

% if null estimates are requested, run stats on all those runs
if ~isempty(nperm)
    
    disp('Estimating the null network statistics...');
    [ nglb, nnde, nnts ] = cellfun(@runStatistics, null, 'UniformOutput', false);
    
    % estimate the pvalues of global measures
    for meas = 1:length(gmeas)
        gbii = gmeas{meas}; % pull the name
        trgb = repmat({glob.(gbii).obs},[nperm 1]); % pull the repeated true value
        
        % estimate and store the p-value
        glob.(gbii).pval = sum(cellfun(@(x,y) x.(gbii).obs >= y, nglb, trgb)) / nperm;
    end
    
    clear meas gbii trgb
    
    % estimate the pvalues of nodal measures
    for meas = 1:length(nmeas)
        noii = nmeas{meas}; % pull the name
        trno = repmat({node.(noii).obs},[nperm 1]); % pull the repeated true values
        
        % pull the null node values greater than truth
        z1 = cellfun(@(x,y) x.(noii).obs >= y, nnde, trno, 'UniformOutput', false);
        
        % estimate and store the node p-values
        node.(noii).pval = sum(cell2mat(z1'), 2) ./ nperm;
    end
    
    % estimate the pvalues of edge measures
    for meas = 1:length(emeas)
        edii = emeas{meas}; % pull the edge wise name
        tred = repmat({nets.(edii).obs},[nperm 1]); % pull the repeated true value
        
        % pull the matrices summarizing the null values
        tre1 = cellfun(@(x,y) x.(edii).obs > y, nnts, tred, 'UniformOutput', false);
        
        % combine null networks into a single matrix
        tre2 = cat(3, tre1{:});
        
        % estimate and store the edge p-values
        nets.(edii).pval = (sum(tre2, 3) >= nets.(edii).obs) ./ nperm;
    end
    
end

%% store statistic estimate back in netw

% store global statistics
netw.stats = glob;

% store nodal statistics in each node
for ii = 1:nnodes
    for jj = 1:size(nmeas, 1)
        netw.nodes{ii}.stats.(nmeas{jj}).obs = node.(nmeas{jj}).obs(ii);
        try
            netw.nodes{ii}.stats.(nmeas{jj}).pval = node.(nmeas{jj}).pval(ii);
        catch
            netw.nodes{ii}.stats.(nmeas{jj}).pval = [];
        end
    end
end

% store edge statistics in each edge
for ii = 1:nedges
    for jj = 1:size(emeas, 1)
        netw.edges{ii}.stats.(emeas{jj}).obs = nets.(emeas{jj}).obs(ii);
        try
            netw.edges{ii}.stats.(emeas{jj}).pval = nets.(emeas{jj}).pval(ii);
        catch
            netw.edges{ii}.stats.(emeas{jj}).pval = [];
        end
    end
end

% store null matrix
if ~isempty(nperm)
    netw.stats.null_mats = cat(3, null{:});
end

end

% this doesn't need to be an internal function, does it?
% as long as globe/node/nets are structs with x.meas.obs/x.meas.pval 

% the function that estimates the stats and return structs of glob/node/nets
function [ glob, node, nets ] = runStatistics(mat)

% a node measure
node.degree.obs = degrees_und(mat)';
node.degree.pval = [];

% a global measure
glob.glbEff.obs = efficiency_wei(mat, 0);
glob.glbEff.pval = [];

% a nets measure
nets.len.obs = weight_conversion(mat, 'lengths');
nets.len.pval = [];

end