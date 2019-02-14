function [ jout ] = fnEdgesToJson(pconn, fg, outdir)
%[ jout ] = fnEdgesToJson(pconn) creates a json file to render the created network.
%
%   Detailed explanation goes here
%
% Brent McPherson, (c) 2019, Indiana University
%

%% combine streamline data from every 50 edges

jj = 1;
count = 1;
ofib = {};

% for every edge
for ii = 1:size(pconn, 1)
    
    % reset count
    if jj > 50
        ofib{count} = coords;
        count = count + 1;
        jj = 1;
        clear coords
    end
    
    % store the nodes indices
    jout(ii).roi1 = pconn{ii}.roi1; % make integer if its not
    jout(ii).roi2 = pconn{ii}.roi2;
    jout(ii).weights.count = pconn{ii}.all.matrix.count;
    jout(ii).weights.density = round(pconn{ii}.all.matrix.density, 3);
    
    % grab the streamlines
    tcoord = fg.fibers(pconn{ii}.all.indices);
    
    if size(tcoord, 1) > 0
        for kk = 1:size(tcoord, 1)
            coords{jj}{kk} = round(tcoord{kk}, 1);
        end
    end
    
    % create the filename
    tname = [ 'conn_' num2str(count) '.json' ];
    
    % if it's empty, no file name
    if isempty(tcoord)
        jout(ii).filename = "";
        jout(ii).idx = "";
    else
        jout(ii).filename = tname;
        jout(ii).idx = jj;    
        jj = jj + 1;
    end
    
end

clear ii jj kk
ofib = ofib';

for ii = 1:size(ofib, 1)
    
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

% replace missing / null w/ zero to fill in empty connections
% make sure 0 is an invalid index or use -1

% save the json out
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
