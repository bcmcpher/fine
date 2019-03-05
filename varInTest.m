function [ out ] = varInTest(parc, fg, names, varargin)
%
% develop the variable input arguments
% 

if mod(size(varargin, 2), 2) == 1
    error('Optional arguments must be passed as named pairs.');
end

% grab name / value pairs
nam = varargin(1:2:end);
val = varargin(2:2:end);

% make sure all names are strings
all(cellfun(@(x) ischar(x), nam));

% make sure all data arrays are numeric and the right size
nfib = size(fg.fibers, 1);
all(cellfun(@(x) size(x, 1) == nfib, val));
all(cellfun(@(x) isnumeric(x), val));

disp([ 'Parcellation file name: ' parc.fname ]);
disp([ 'Number of fg streamlines: ' num2str(size(fg.fibers, 1)) ]);
disp([ 'Names vector size: ' num2str(size(names, 1)) ]);

for ii = 1:size(nam, 2)
    out.(nam{ii}) = val{ii};
end

end

