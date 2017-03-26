function [ hmap ] = plot_endpoint_heatmap(fe, outfile)
%% end point density plot
%
% creates a heatmap of the streamline terminations in DWI image space
% matches ACPC alignment, in register w/ anat
%
% add ability to get specific ROIs?
% create reliced count in anat space?
%

% initialize endpoint outputs
iep1 = zeros(length(fe.fg.fibers), 3);
iep2 = zeros(length(fe.fg.fibers), 3);

% for every fiber, pull the end points
for ii = 1:length(fe.fg.fibers)
    iep1(ii,:) = fe.fg.fibers{ii}(:,1)';
    iep2(ii,:) = fe.fg.fibers{ii}(:,end)';
end

% round all the endpoints
iep1 = round(iep1) + 1;
iep2 = round(iep2) + 1;

% combine and reconver fiber endpoints to image space
ep = [ iep1; iep2 ];

% sanity check of low/high dims
%minmax(ep')

% create empty data space
img = zeros(fe.life.imagedim(1:3));

% for every endpoint, create a density map of rois
for ii = 1:length(ep)
    img(ep(ii,1), ep(ii, 2), ep(ii, 3)) = img(ep(ii,1), ep(ii, 2), ep(ii, 3)) + 1;
end

% create output nifti
hmap = niftiCreate('data', img, 'fname', outfile, ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);
% any additional parameters?
               
% save file
niftiWrite(hmap, outfile);

end

