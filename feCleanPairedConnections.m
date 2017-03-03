function [ pconn, cln ] = feCleanPairedConnections(fe, pconn, label)
%feCleanPairedConnections cleans the outlier fibers from a pconn cell
% array. 
%
% INPUTS:
% - 'fe' is the fit fe structure that generated the paired connections object
% - 'pconn' is the paired connections object that has the field to be cleaned
% - 'label' is the is the field of indices that will be cleaned.
%   
% OUTPUTS:
% - 'pconn' is the paired connections object wiht the cleaned streamline
%   indices added in a new field appended with '_clean'.
%
% - 'cln' is the cell array of length(pconn) that is created in the fxn and
%   added to pconn. It is for debugging - no different information is
%   stored here.
%

display('Converting streamlines to ACPC space...');

% create fiber group
fg = feGet(fe, 'fg acpc');

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

display(['Started cleaning ' num2str(length(pconn)) ' paired connections...']);

tic;
parfor ii = 1:length(pconn)
    
    % get the requested field
    tmp = getfield(pconn{ii}, label);
    
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
        
        % create an fg group of the connection
        cln{ii}.connfib = fgExtract(fg, cln{ii}.indices, 'keep');
        
        % for all streamlines, compute outliers
        [ ~, cln{ii}.all ] = mbaComputeFibersOutliers(cln{ii}.connfib, maxDist, maxLengthStd, numNodes, 'mean', 0, maxIter);

        % catch outputs
        cln{ii}.out.indices = cln{ii}.indices(cln{ii}.all);
        cln{ii}.out.lengths = cln{ii}.lengths(cln{ii}.all);
        cln{ii}.out.weights = cln{ii}.weights(cln{ii}.all);
                
        % catch image coords for link network
        %cln{ii}.croi = fefgGet(cln{ii}.afg, 'unique image coords');

%         % I don't clip fibers or split loops - but do I really need to?
%         % I've already got the end point indices, do I have to clean this way?
%         % That approach is more targeted AFQ, not end point operations
%        
%         % create connection fg structure
%         cln{ii}.connfib = fgExtract(fg, cln{ii}.indices, 'keep');
%         
%         % pull path neighborhood connection and ROI
%         cln{ii}.pnfibs = feGet(fe, 'Path Neighborhood', cln{ii}.indices);
%         cln{ii}.indroi = feGet(fe, 'coords from fibers', cln{ii}.pnfibs);
%         
%         % initial cleaning
%         cln{ii}.clip = feClipFibersToVolume(cln{ii}.pnfibs, cln{ii}.indroi, maxVolDist);
%         cln{ii}.clip.fibers = mbaFiberSplitLoops(cln{ii}.clip.fibers);
%         
%         % remove fibers below a certain length
%         % how can I keep indices tracked through here?
%         
%         % compute outliers
%         [ ~, keep ] = mbaComputeFibersOutliers(tafg, 3, 3, 100, 'mean', 0, 10);
     
        % keep count of cleaned connections
        clncnt = clncnt + 1;
        clntot(ii) = sum(cln{ii}.all);
        
    else 
        % fill in empty cells and skip to next connection
        cln{ii}.out.indices = [];
        cln{ii}.out.lengths = [];
        cln{ii}.out.weights = [];
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
for ii = 1:length(pconn)

    % values used to calculate newly cleaned connections
    psz = pconn{ii}.roi1sz + pconn{ii}.roi2sz;
    cnt = size(cln{ii}.out.indices, 1);
    len = mean(cln{ii}.out.lengths);
    if isempty(cln{ii}.out.lengths) % if there are no nz lengths
        dln = 0;
    else
        dln = sum(1 / cln{ii}.out.lengths);
    end
    
    % create all cleaned streamline counts
    matrix.count = cnt;
    matrix.density = (2 * cnt) / psz;
    matrix.length = len;
    matrix.denlen = (2 / psz) * dln;
    
    % create structure with all the cleaned field data
    tmp = struct('indices', cln{ii}.out.indices, ...
                 'lengths', cln{ii}.out.lengths, ...
                 'weights', cln{ii}.out.weights, ...
                 'matrix', matrix);
        
    % assign cleaned count matrix
    pconn{ii} = setfield(pconn{ii}, newout, tmp);
    
end
    
end
