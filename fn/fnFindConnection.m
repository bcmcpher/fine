function [ conn, out ] = fnFindConnection(fg, pconn, label, roi1, roi2)
%fnFindConnection() returns a fiber group between to desired indices if one
% exists.
%
% add profile, matrix, etc. outputs
%

% sort ROI indices
rois = sort([ roi1 roi2 ]);

% for every connection
for conn = 1:size(pconn, 1)
    
    % find the first and second
    if pconn{conn}.roi1 == rois(1)
        if pconn{conn}.roi2 == rois(2)
            
            % grab the labels indices
            ind = pconn{conn}.(label).indices;
            
            if isempty(ind)
                error('The requested connection is empty.');
            end
            
            % create the output connection and length
            out = fgCreate('name', [ 'fg_' num2str(roi1) '_' num2str(roi2) ], ...
                           'colorRgb', [ .1 .25 .65 ], 'fibers', fg.fibers(ind));
            length = pconn{conn}.(label).matrix.length;
                       
            % leave the loop
            break
        end
    end

end

% print some information
display(['Connection between ROI''s numbered ' num2str(roi1) ' and ' num2str(roi2) ' found at index: ' num2str(conn) ]);
display(['The fascicle has ' num2str(size(ind, 1)) ' fasciles and an average length of ' num2str(length) ' mm.']);

end