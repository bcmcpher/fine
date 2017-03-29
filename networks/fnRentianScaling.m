function [ p, se, fh ] = fnRentianScaling(mat, rois, reps)
%fnRentianScaling compute and display network measure that controls for
% the physical space of the network.
% 

% binarize input
dat = weight_conversion(mat, 'binarize');

% grab node coords from rois object 
coords = zeros(length(rois), 3);
for ii = 1:length(rois)
    coords(ii, :) = rois{ii}.centroid.acpc;
end

% rentian scaling - separate fxn
[ N, E ] = rentian_scaling_3d(dat, coords, reps, 0.000001); % n is number of partitions

% estimate Rent's exponent
[ b, stats ] = robustfit(log10(N), log10(E));

% return Rent's exponent and standard error from robustfit()
% THIS IS COMPLETELY DIFFERENT FROM DESCRIPTION, MY BEST GUESS
p = b(2);
se = stats.se(2);

% plot from description...
fh = figure; 
loglog(E, N, '*');

end

