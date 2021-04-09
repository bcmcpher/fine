function [ netw ] = fnCleanEdges(netw, fg, minLength, maxDist, maxLengthStd, numNodes, maxIter, minNum)
%feCleanPairedConnections cleans the outlier fibers from connections inside 
% a paired connection (pconn) cell array. 
%
% INPUTS:
%    fg        - fiber group in acpc space
%    pconn     - paired connection object created with fg
%    label     - string indicating the fiber groups to clean
%                either:
%                       'all' for all assigned streamlines or
%                       'nzw' for non-zero weighted fibers returned by LiFE
%    minLength    - the minimum length in mm of a connection to be considered
%                   valid (default = 10)
%    maxDist      - the # of standard deviations a streamline node can take from the
%                   centroid of the tract (default = 4)
%    maxLengthStd - the # of standard deviations the length of a
%                   streamline can deviate from the tracts mean (default = 4)
%    numNodes     - the number of nodes fibers are resampled to for
%                   similarity comparisons (default = 100)
%    maxIter      - the number of iterations outliers will be tested and
%                   sampled from (default = 10)
%    minNum       - the minimum number of streamlines that must exist for
%                   a connection to be kept (default = 3).
%                   NOTE: tract profiles require a minimum of 3 streamlines.
%             
% OUTPUTS:
%     pconn - is the paired connections object with the cleaned streamline
%             indices added in a new label appended with '_clean'.
%
% EXAMPLE:
%
% % load data
% parc = niftiRead('labels.nii.gz');
% fg        = feGet(fe, 'fibers acpc');
% fibers    = fg.fibers;
% fibLength = fefgGet(fg, 'length');
% weights   = feGet(fe, 'fiberweights');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % clean streamlines too far from center with all streamlines that are assigned
% pconn = feCleanPairedConnections(fg, pconn, 'all');
%
% Brent McPherson (c), 2017 - Indiana University
%

% fill in defaults for any missing parameters

if(~exist('minLength', 'var') || isempty(minLength))
    minLength = 10;
end

if(~exist('maxDist', 'var') || isempty(maxDist))
    maxDist = 4;
end

if(~exist('maxLengthStd', 'var') || isempty(maxLengthStd))
    maxLengthStd = 4;
end

if(~exist('numNodes', 'var') || isempty(numNodes))
    numNodes = 100;
end

if(~exist('maxIter', 'var') || isempty(maxIter))
    maxIter = 10;
end

if(~exist('minNum', 'var') || isempty(minNum))
    minNum = 3;
end

% add cleaning parameters to output
netw.clean.minLength = minLength;
netw.clean.maxDist = maxDist;
netw.clean.maxLengthStd = maxLengthStd;
netw.clean.numNodes = numNodes;
netw.clean.maxIter = maxIter;
netw.clean.minNum = minNum;

% pull a logical from each edge if any streamlines exist
pc = sum(cellfun(@(x) size(x.fibers.indices, 1) > 0, netw.edges));
disp(['Started cleaning ' num2str(pc) ' present connections...']);

% initialize a count of before / after cleaning edges
bfcln = zeros(size(netw.edges, 1), 1);
afcln = zeros(size(netw.edges, 1), 1);

% pull connections from network
edges = netw.edges;

tic;
for ii = 1:size(edges, 1)
    
    % get the requested edge
    edge = edges{ii};
    
    % grab the before count
    bfcln(ii) = size(edge.fibers.indices, 1);    
    
    % if the connection isn't empty
    if ~isempty(edge.fibers.indices)
                
        % drop indices less than a particular length
        keep_idx = edge.fibers.lengths > minLength;
                
        % if there are too few streamlines, fill in empty and move on
        if(sum(keep_idx) < minNum)
            edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
            edges{ii} = edge;
            continue
        else
           % filter other values before fg cleaning
           edge.fibers = structfun(@(x) x(keep_idx), edge.fibers, 'UniformOutput', false);
        end
        
        % create an fg group of the filtered edge
        tfg = fgExtract(fg, edge.fibers.indices, 'keep');
        
        % try and clean the pruned length / count, fail over to empty 
        try
            % for all streamlines, compute outliers
            % cln is a logical of the streamlines to keep
            [ ~, cln ] = mbaComputeFibersOutliers(tfg, maxDist, maxLengthStd, numNodes, 'mean', 0, maxIter);
            
            % filter edge down to cleaned fg
            edge.fibers = structfun(@(x) x(cln), edge.fibers, 'UniformOutput', false);
            
        catch
            warning('Cleaning of edge index %d failed. Assigning empty values in network structure.', ii);
            edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
            edges{ii} = edge;
            continue
        end
        
        % keep count of cleaned connections
        % (before cleaning) - (total > minLength and after cleaning fxn)
        afcln(ii) = bfcln(ii) - sum(cln);
        
    else 
        
        % fill in empty cells and skip to next connection
        edge.fibers = structfun(@(x) [], edge.fibers, 'UniformOutput', false);
        edges{ii} = edge;
        continue
        
    end
    
    % just in case it's missed, just fill it back in
    edges{ii} = edge;
    
end
time = toc;

clear ii edge keep_idx tfg cln

% reassign paired connection
netw.edges = edges;

% compute how many are kept vs. dropped
diffs = bfcln - afcln;

disp(['Cleaned outlier streamlines from ' num2str(sum(afcln ~= bfcln)) ' edges in ' num2str(round(time)/60) ' minutes.']);
disp(['Removed ' num2str(sum(afcln)) ', keeping ' num2str(sum(diffs)) ' streamlines of ' num2str(sum(bfcln)) '.']);

% track total streamlines dropped
netw.clean.drop = sum(afcln);
netw.clean.diffs = afcln;
    
end
