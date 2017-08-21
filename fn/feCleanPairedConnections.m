function [ pconn, cln ] = feCleanPairedConnections(fg, pconn, label)
%feCleanPairedConnections cleans the outlier fibers from connections inside 
% a paired connection (pconn) cell array. 
%
% INPUTS:
%     fg    - fiber group in acpc space
%     pconn - paired connection object created with fg
%     label - string indicating the fiber groups to clean
%             either:
%                     'all' for all assigned streamlines or
%                     'nzw' for non-zero weighted fibers returned by LiFE
%             
% OUTPUTS:
%     pconn - is the paired connections object with the cleaned streamline
%             indices added in a new label appended with '_clean'.
%
%     cln   - debugging output; the cell array that is added internally to pconn
%
% TODO:
% - provide fxn level access to cleaning parameters
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

% preallocate cleaned edge structure
cln = cell(length(pconn), 1);
clncnt = 0;
clntot = 0;

% cleaning arguments - set up as arguments
%maxVolDist = 3;
minLength = 10;
maxDist = 4;
maxLengthStd = 4;
numNodes = 100;
maxIter = 10;
minNum = 3;

display(['Started cleaning ' num2str(length(pconn)) ' paired connections...']);

tic;
parfor ii = 1:length(pconn)
    
    % get the requested field
    tmp = pconn{ii}.(label);
    
    if ~isempty(tmp.indices)
        
        % pull indices of connection
        cln{ii}.indices = tmp.indices;
        cln{ii}.lengths = tmp.lengths;
        cln{ii}.weights = tmp.weights;
        
        % drop indices less than a particular length
        cln{ii}.lenindx = cln{ii}.lengths > minLength;
        cln{ii}.indices = cln{ii}.indices(cln{ii}.lenindx);
        cln{ii}.lengths = cln{ii}.lengths(cln{ii}.lenindx);
        cln{ii}.weights = cln{ii}.weights(cln{ii}.lenindx);
        
        % if there are too few streamlines, fill in empty and move on
        if(size(cln{ii}.indices, 1) < minNum)
            cln{ii}.out.indices = [];
            cln{ii}.out.lengths = [];
            cln{ii}.out.weights = [];
            %cln{ii}.out.pvoxels = [];
            continue
        end
        
        % create an fg group of the connection
        cln{ii}.connfib = fgExtract(fg, cln{ii}.indices, 'keep');
        
        % for all streamlines, compute outliers
        [ ~, cln{ii}.all ] = mbaComputeFibersOutliers(cln{ii}.connfib, maxDist, maxLengthStd, numNodes, 'mean', 0, maxIter);

        % catch outputs
        cln{ii}.out.indices = cln{ii}.indices(cln{ii}.all);
        cln{ii}.out.lengths = cln{ii}.lengths(cln{ii}.all);
        cln{ii}.out.weights = cln{ii}.weights(cln{ii}.all);
        
        % keep count of cleaned connections
        clncnt = clncnt + 1;
        clntot(ii) = sum(cln{ii}.all);
        
    else 
        
        % fill in empty cells and skip to next connection
        cln{ii}.out.indices = [];
        cln{ii}.out.lengths = [];
        cln{ii}.out.weights = [];
        %cln{ii}.out.pvoxels = [];
        continue
        
    end
    
end
time = toc;

display(['Cleaned outlier streamlines from ' num2str(clncnt) ' edges in ' num2str(round(time)/60) ' minutes.']);
display(['Kept ' num2str(sum(clntot)) ' streamlines.']);
clear ii tmp

% create new output in pconn
newout = strcat(label, '_clean');

% add cleaned streamlines to pconn
parfor ii = 1:length(pconn)

    % values used to calculate newly cleaned connections
    psz = pconn{ii}.roi1sz + pconn{ii}.roi2sz;
    cnt = size(cln{ii}.out.indices, 1);
    len = mean(cln{ii}.out.lengths);
    if isempty(cln{ii}.out.lengths) % if there are no nz lengths
        dln = 0;
    else
        dln = sum(1 / cln{ii}.out.lengths);
    end
    
    % create structure with all the cleaned field data
    tmp = struct('indices', cln{ii}.out.indices, ...
                 'lengths', cln{ii}.out.lengths, ...
                 'weights', cln{ii}.out.weights);
             
    % assign cleaned count matrix
    pconn{ii}.(newout) = tmp;
    pconn{ii}.(newout).matrix.count = cnt;
    pconn{ii}.(newout).matrix.density = (2 * cnt) / psz;
    pconn{ii}.(newout).matrix.length = len;
    pconn{ii}.(newout).matrix.denlen = (2 / psz) * dln;
    
end
    
end
