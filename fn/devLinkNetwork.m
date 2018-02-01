function [ omat, time ] = devLinkNetwork(pconn, label, Phi, dict)
%devLinkNetwork Summary of this function goes here
%   Detailed explanation goes here

%% generate link network

% error if volume is not precomputed
if (~isfield(pconn{1}.(label), 'volume'))
    error('Edge volumes for this label must be precomputed with fnFindPathVoxels().');
end    

disp('Preallocating output matrix...');

% number of edges; links between edges are the new nodes
lnodes = size(pconn, 1);

% compute indices of output size
mask = triu(ones(lnodes), 1); 

% preallocate unique indices - faster than nchoosek
[ xind, yind ] = find(mask > 0);
indx = [ xind, yind ];

clear mask xind yind 

% build an empty output array / matrix
out = nan(size(indx, 1), 1);
omat = nan(lnodes, lnodes);

disp('Computing weights for unique edge intersections...');

tic;
parfor ii = 1:size(indx, 1)
    
    % grab the precomputed indices
    conn1 = indx(ii, 1);
    conn2 = indx(ii, 2);

    % grap the unique streamlines of each edge
    pi1 = pconn{conn1}.(label).indices;
    pi2 = pconn{conn2}.(label).indices;
    
    % fill in nothing if either connection is empty
    if isempty(pi1) || isempty(pi2)
        out(ii) = 0;
        continue
    end
    
    % grab the unique voxels of each edge
    pv1 = pconn{conn1}.(label).volume.pvoxels;
    pv2 = pconn{conn2}.(label).volume.pvoxels;
        
    % find the common volume of the connections
    vind = intersect(pv1, pv2);
    
    % if the intersection is empty, fill in zero
    if isempty(vind)
        out(ii) = 0;
        continue
    end
    
%     % dice coefficient
%     
%     % find the size of the intersection - numerator
%     num = size(vind, 1);
%     
%     % combine the voxel indices for denominator of Dice coeff
%     den = size(pv1, 1) + size(pv2, 1);
%     
%     % build Dice coeff values for assignment
%     out(ii) = (2 * num) / den;
    
    % angle of intersection (subset the tensor)
    
    % reduce Phi to the shared voxels of the connections
    sub1 = Phi(:, vind, pi1);
    sub2 = Phi(:, vind, pi2);
    
    % from subtensors pull the dictionary atoms for each connection
    atom1 = dict(:, sub1.subs(:, 1));
    atom2 = dict(:, sub2.subs(:, 1));
    
    % find the average orientation of each tract
    matm1 = mean(atom1, 2);
    matm2 = mean(atom2, 2);
    
    % find the scalar product to assign a link weight
    out(ii) = dot(matm1, matm2);
    
    % other tensor values that could be returned
    
    % compute the angle between connections
    %theta = acos(val);
    
    % add switch to either compute mean or use pdist2 on atom values
    %csim = pdist2(atom1, atom2, 'cosine'); % is not a single value?
    % how did I use this - or something - to get a single value?
    
end
time = toc;

disp(['Computed link network of intersecting angles in ' num2str(time/60) ' minutes.' ]);

disp('Creating link network...');

% create outputs
for ii = 1:size(indx, 1)
    
    % grab the precomputed indices
    conn1 = indx(ii, 1);
    conn2 = indx(ii, 2);
    
    % create matrix
    omat(conn1, conn2) = out(ii);
    omat(conn2, conn1) = out(ii);
    
end

% fix nan/inf/neg values to zero
omat(isinf(omat)) = 0;
omat(isnan(omat)) = 0;
%omat(omat < 0) = 0;

end

