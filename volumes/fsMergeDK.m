function [ names, olabs ] = fsMergeDK(aparc, outfile)
%% 
% currently assigns mode value of neighbors
% can repeat inflation, unsure of ideal optimum number of times
% check to see if I can import FS directly?
%
% aparc = '105115/mri/aparc+aseg.nii.gz'; infl = 3;

%% hardcoded labels for this helper conversion

% predefined cortex / subcortical labels for this parcellation
ctxvals = double([ 1001:1035 2001:2035 ])';
sbcvals = [ 10 11 12 13 17 18 26 49 50 51 52 53 54 58 ]';

% hard code labels into lobes for a courser parcellation
lfrontal = [ 1003 1012 1014 1017 1018 1019 1020 1024 1027 1028 1032 ];
lparietal = [ 1008 1022 1025 1029 1031 ];
ltemporal = [ 1001 1006 1007 1009 1015 1016 1030 1033 1034 1035 ];
loccipital = [ 1005 1011 1013 1021 ]; 
lcingulate = [ 1002 1010 1023 1026 ];

% increase labels for right labels
rfrontal = lfrontal + 1000;
rparietal = lparietal + 1000;
rtemporal = ltemporal + 1000;
roccipital = loccipital + 1000;
rcingulate = lcingulate + 1000;

% merge into a list so the loop can replace the indices in each entry
lobes = { lfrontal, lparietal, ltemporal, loccipital, lcingulate, rfrontal, rparietal, rtemporal, roccipital, rcingulate };
names = {'Left Frontal', 'Left Parietal', 'Left Temporal', 'Left Occipital', 'Left Cingulate', ...
         'Right Frontal', 'Right Parietal', 'Right Temporal', 'Right Occipital', 'Right Cingulate'};
     
%% process input
disp('Loading aparc+aseg nifti...');

% read in aligned aparc+aseg .nii (read labels)
rlabs = niftiRead(aparc);
rdat = rlabs.data;

% pull the labels from the input volume to check that they exist
ctxrval = unique(rlabs.data(:));

% check the size of the input labels present in the volume
if size(intersect(ctxrval, ctxvals), 1) ~= 68
    warning('Some expected labels are not present in volume.\nThis fxn will likely not perform as expected.');
end

% check if any subcortical labels are present
% if they are, add the subcortical labels to the merged parcellation
if size(intersect(ctxrval, sbcvals), 1) > 0
    lobes = [ lobes sbcvals ];
    names = [ names 'Subcortical' ];
end

% start relabeling volume at 901
fix = 901;

% for every major lobe
for lobe = 1:length(lobes)
    
    % pull all the labels to rename
    labs = lobes{lobe};
    
    % for every label in the lobe
    for lab = 1:length(labs)
        
        rdat(rdat == labs(lab)) = fix;
        
    end
    
    % iterate new labels
    fix = fix + 1;

end
    
% create merged labels nifti
olabs = rlabs;
olabs.fname = outfile;
olabs.data = rdat;

disp(['Saving merged ROIs in: ' olabs.fname '...']);

% save inflated labels
niftiWrite(olabs, olabs.fname);

end
