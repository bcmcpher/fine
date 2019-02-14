function [ hmap ] = fnEndpointFg(fg, anat, outfile, smooth)
%hmap = fnEndpointHeatmap() creates a heatmap on an fg structures endpoints.
% sliced to a specified volume and optionally smoothed. 
% It returns a nifti object with the extracted data.
% 
% INPUTS:
%     fg      - a fiber group in AC-PC space
%     anat    - the volume space to export data into
%     outfile - the name of the .nii.gz file to write the volume to
%     smooth  - perform 3D smoothing of provided kernel size
% 
% OUTPUTS:
%     hmap    - a .nii structure containing the computed heatmap data
%
% Brent McPherson (c) 2019 Indiana University
%

%% parse optional arguments

disp(['Image will be resliced into ' anat.fname ' space...']);

out_xyz = anat.qto_xyz;
out_ijk = anat.qto_ijk;
out_dim = size(anat.data);    

% turn off smoothing if it's not requested
if(~exist('smooth', 'var') || isempty(smooth))
    smooth = [];
    disp('Not smoothing image...');
else
    disp(['Image will be smoothed with a [ ' num2str(smooth) ' ] gaussian kernel...']);
end

disp('Converting streamlines to output space...');

% convert to output space
fg = dtiXformFiberCoords(fg, out_ijk, 'img');

% initialize endpoint outputs
iep1 = zeros(length(fg.fibers), 3);
iep2 = zeros(length(fg.fibers), 3);

% for every fiber, pull the end points
for ii = 1:length(fg.fibers)
    iep1(ii,:) = fg.fibers{ii}(:,1)';
    iep2(ii,:) = fg.fibers{ii}(:,end)';
end

% combine fiber endpoints in acpc space
ep = [ iep1; iep2 ];
ep = round(ep) + 1;

clear iep1 iep2

disp('Creating output data...');

% create empty data space
img = zeros(out_dim);

% create a density map of rois from every termination
for ii = 1:length(ep)
    img(ep(ii,1), ep(ii, 2), ep(ii, 3)) = img(ep(ii, 1), ep(ii, 2), ep(ii, 3)) + 1;
end

% smooth data if requested
if ~isempty(smooth)
    disp('Smoothing data...');
    img = smooth3(img, 'gaussian', smooth);
end

% create output nifti - write to fe space or new space
hmap = niftiCreate('data', img, 'fname', outfile, ...
                   'qto_xyz', out_xyz, ...
                   'qto_ijk', out_ijk);
               
% save file
niftiWrite(hmap, outfile);

end

