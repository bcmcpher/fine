function [ hmap ] = plot_endpoint_heatmap(fe, outfile, smooth, reslice)
%% end point density plot
%
% creates a heatmap of the streamline terminations from an fe structure.
% The data can optionally be smoothed with a 3d gaussian or resliced into a
% different image space. Reslicing to a new space should proabably be smoothed.
%
%

% parse optional arguments

% turn off smoothing if it's not requested
if(~exist('smooth', 'var') || isempty(smooth))
    smooth = [];
    disp('Not smoothing image...');
end

% if no new space is provided, use the fe structure for output size
if(~exist('reslice', 'var') || isempty(reslice))
    disp('Image will be sliced in fe space...');
    out_xyz = fe.life.xform.img2acpc;
    out_ijk = fe.life.xform.acpc2img;
    out_dim = fe.life.imagedim(1:3);

% otherwise use the data from the new space for output size
else
    disp(['Image will be resliced into ' reslice.fname ' space...']);
    out_xyz = reslice.qto_xyz;
    out_ijk = reslice.qto_ijk;
    out_dim = size(reslice.data);    
end

disp('Converting streamlines to AC-PC space and extracting end points...');

% pull acpc fibers
fg = feGet(fe, 'fibersacpc');

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
    disp(['Smoothing output with a [ ' num2str(smooth) ' ] gaussian kernel...']);
    img = smooth3(img, 'gaussian', smooth);
end

% create output nifti - write to fe space or new space
hmap = niftiCreate('data', img, 'fname', outfile, ...
                   'qto_xyz', out_xyz, ...
                   'qto_ijk', out_ijk);
               
% save file
niftiWrite(hmap, outfile);

end

