function [ LNmap, MLmap ] = fnLinkNetworkVolumes(fe, pconn, label, assign, one_out, mlt_out, minVol)
%fnLinkNetworkVolumes converts link network communities to voxel volumes
%   with some additional work, I hope this will identify major tracts
%
% [ LNmap, MLmap ] = fnLinkNetworkVolumes(fe, pconn, 'nzw', node.assign, 'zz_one.nii.gz', 'zz_two.nii.gz');
%

% define defaults for minimum common voxels
if(~exist('minVol', 'var') || isempty(minVol))
    minVol = 1;
end

% grab count of each community
[ cnt, lab ] = hist(assign, unique(assign));
cnt = cnt';

% create structure of label and size of each link community
linkLabelSize = [ lab, cnt ];

% show the linkLabels that are greater than 1
display(['Found: ' num2str(size(linkLabelSize(cnt > 1, :), 1)) ' labels.']);

% find the actual connections worth making volumes of
linkLabels = linkLabelSize(cnt > 1, 1);

% preallocate output
out = cell(size(linkLabels, 1), 1);

% loop over labels, grab voxel indices into tensor
for ii = 1:size(linkLabels, 1)
    
    % grab label here
    lab = linkLabels(ii);
    
    % create entry
    out{ii}.label = lab;
    out{ii}.voxel = [];
    
    % for every connection
    for jj = 1:size(pconn)
        
        % if the edge assignment is the label of interest
        if assign(jj) == lab
            
            % append to voxels to that ROI
            out{ii}.voxel = [ out{ii}.voxel; pconn{jj}.(label).volume.pvoxels ];
            
        end
        
    end
    
    % only grab voxels that occur more than once
    vol = unique(out{ii}.voxel);
    num = numel(vol);
    count = zeros(num,1);
    
    % for every unique voxel
    for k = 1:num
        
        % get the count
        count(k) = sum(out{ii}.voxel == vol(k));
        
    end
    
    % debug printing max intersection of communities
    display([ 'Label ' num2str(ii) ' has as many as ' num2str(max(count)) ' shared voxels.' ]);
    
    % set the unique voxels across label
    out{ii}.voxel = vol(count > minVol);
    %out{ii}.voxel = unique(out{ii}.voxel);
    
end

% create output structure to loop over and fill in volumes

% create image size
imgsize = fe.life.imagedim(1:3);

% create empty image of zeros that matches coordinate space
oneOut = zeros([ imgsize 1 ]);
mltOut = zeros([ imgsize size(out, 1) ]);

% pull fe roi volume
roi = feGet(fe, 'roicoords');

% for every label
for ii = 1:size(out, 1)
    
    % pull coordinates from ROI index into volume
    coord = roi(out{ii}.voxel, :);

    % for every unique voxel in the link community
    for jj = 1:size(coord, 1)
    
        % write out to a single, multi-label volume
        oneOut(coord(jj, 1), coord(jj, 2), coord(jj, 3)) = ii;
        
        % write out to a multi-volume mask
        mltOut(coord(jj, 1), coord(jj, 2), coord(jj, 3), ii) = ii;
    
    end
    
end

% create output niftis
LNmap = niftiCreate('data', oneOut, 'fname', one_out, ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);

MLmap = niftiCreate('data', mltOut, 'fname', mlt_out, ...
                   'qto_xyz', fe.life.xform.img2acpc, ...
                   'qto_ijk', fe.life.xform.acpc2img);

% save down niftis
niftiWrite(LNmap, one_out);
niftiWrite(MLmap, mlt_out);

end

