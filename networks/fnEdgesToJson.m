function [ jout, ofib ] = fnEdgesToJson(pconn, fg, outdir)
%[ jout, ofib ] = fnEdgesToJson(pconn, fg, outdir) creates a json file to render the created network.
%
%   This function creates the output folder of .json files needed for the
%   visualization brainlife.io to render the network edge weights alongside
%   the anatomy.
%
% INPUTS:
%     pconn  - paired connection object that stores network edge weights
%              and streamline indices
%     fg     - fiber group that corresponds to streamline indices in pconn
%     outdir - the output directory where the .json files will be written
%
% OUTPUTS:
%     jout - the index.json file that stores edge properties
%     ofib - stores the streamline data for export to files
%            each cell has streamline data for 50 edges
%
% This is still being actively developed, some features may change.
%
% Brent McPherson, (c) 2019, Indiana University
%

% initialize counter to store streamline data by groups of 50
jj = 1;
count = 1;
ofib = {};

disp('Indexing network data and streamlines...');

% for every edge
for ii = 1:size(pconn, 1)
    
    % once 50 have been stored
    if jj > 50
        % iterate the object / reset the count
        count = count + 1;
        jj = 1;
        clear coords
    end
    
    % store the nodes indices
    jout(ii).roi1 = pconn{ii}.roi1;
    jout(ii).roi2 = pconn{ii}.roi2;
    
    % grab every edge weight that is found (need to loop / generalize)
    jout(ii).weights.count = pconn{ii}.all.matrix.count;
    jout(ii).weights.density = round(pconn{ii}.all.matrix.density, 5);
    jout(ii).weights.length = round(pconn{ii}.all.matrix.length, 5);
    jout(ii).weights.denlen = round(pconn{ii}.all.matrix.denlen, 5);
    
    % grab the streamlines from fg
    tcoord = fg.fibers(pconn{ii}.all.indices);
    
    % if the connection exists
    if size(tcoord, 1) > 0
        % for every streamline
        for kk = 1:size(tcoord, 1)
            % store the streamline nodes in a nested array
            coords{jj}{kk} = round(tcoord{kk}, 2);
        end
    end
    
    % create the filename to store the streamline data in
    tname = [ 'conn_' num2str(count) '.json' ];
    
    % if it's empty, no file name
    if isempty(tcoord)
        jout(ii).filename = "";
        jout(ii).idx = "";
    else
    % otherwise assign a name and index for the app
        jout(ii).filename = tname;
        jout(ii).idx = jj-1; % already off by one, correct for 2?    
        jj = jj + 1;
    end
    
    % add streamlines to output structure if they exist
    if size(tcoord, 1) > 0
        ofib{count} = coords;
    end
    
end

clear ii jj kk
ofib = ofib';

disp('Writing streamline data to json...');

% for every collection of 50 files
for ii = 1:size(ofib, 1)
    
    % create an output file name
    tname = [ 'conn_' num2str(ii) '.json' ];
    
    % write it out
    opt.filename = strcat(outdir, '/', tname);
    opt.ArrayIndent = 1;
    opt.ArrayToStruct = 0;
    opt.SingleArray = 0;
    opt.SingletCell = 0;
    opt.Compact = 1;
    savejson('', ofib{ii}, opt);

end

disp('Writing index.json file...');

% save the index (edge weights, etc.) json out
opt.filename = strcat(outdir, "/index.json");
opt.ArrayIndent = 1;
opt.ArrayToStruct = 0;
opt.SingleArray = 0;
opt.SingletCell = 0;
opt.Compact = 1;
savejson('roi_pairs', jout, opt);

end

%x1 = x(:, 1:end-1);
%x2 = x(:, 2:end);
%sum(sqrt(sum((x1 - x2) .^ 2)))
%z = cellfun(@(x) sum(sqrt(sum((x(:, 1:end-1) - x(:, 2:end)) .^ 2))), fg.fibers, 'UniformOutput', true);
