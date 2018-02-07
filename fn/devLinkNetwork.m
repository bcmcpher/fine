function [ omat, time ] = devLinkNetwork(pconn, label, meas, dtype, Phi, dict)
%devLinkNetwork Summary of this function goes here
%   Detailed explanation goes here
%
% TODO:
% - documentation
% - argument parsing / smart outputs
%
%

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
    
    switch meas
        
        % dice coefficient
        case 'dice'
        
            % find the size of the intersection - numerator
            num = size(vind, 1);
            
            % combine the voxel indices for denominator of Dice coeff
            den = size(pv1, 1) + size(pv2, 1);
            
            % build Dice coeff values for assignment
            out(ii) = (2 * num) / den;
        
        case {'mi', 'mutualinfo'}
            
            % SOME INDIVIDUAL ENTROPY ESTIMATES ARE LARGE NEGATIVE NUMBERS
            % THIS MAKES ~HALF OF MI ESTIMATES LARGE NEGATIVE NUMBERS - BAD
            % I DO NOT CURRENTLY KNOW HOW TO FIX THIS
            % link to how I compute the values currently
            % https://stackoverflow.com/questions/23691398/mutual-information-and-joint-entropy-of-two-images-matlab
            
            % grab path voxels and values for the first connection
            imp1 = pconn{conn1}.(label).volume.pvoxels;
            imv1 = pconn{conn1}.(label).volume.(dtype).raw;
            
            % grab path voxels and values for the second connection
            imp2 = pconn{conn2}.(label).volume.pvoxels;
            imv2 = pconn{conn2}.(label).volume.(dtype).raw;
            
            % grab combined voxels from both paths
            imp = union(imp1, imp2);
            
            % build image vector
            im = nan(size(imp, 1), 2);
            
            % grab indices of voxels not in the path
            [ ~, im01 ] = setdiff(imp, imp1);
            [ ~, im02 ] = setdiff(imp, imp2);
            
            % set entries unique to each tract to zero in the other tract
            im(im01, 1) = 0;
            im(im02, 2) = 0;
            
            % replace missing in each column with the values
            im(isnan(im(:, 1)), 1) = imv1;
            im(isnan(im(:, 2)), 2) = imv2;

            % compute bins of joint histogram
            [ ~, ~, indrow ] = unique(im(:, 1));
            [ ~, ~, indcol ] = unique(im(:, 2));
            
            % compute joint entropy
            jointHist = accumarray([indrow indcol], 1);
            jointProb = jointHist / numel(indrow);
            indNoZero = jointHist ~= 0;
            jointNzPb = jointProb(indNoZero);
            jointEntropy = -sum(jointNzPb .* log2(jointNzPb));
            
            % compute individual histogram summaries
            histImage1 = sum(jointHist, 1);
            histImage2 = sum(jointHist, 2);
            
            % find non-zero elements for first image's histogram
            % extract them out and get the probabilities
            % compute the entropy
            indNoZero1 = histImage1 ~= 0;
            prob1NoZero = histImage1(indNoZero1) / numel(histImage1);
            entropy1 = -sum(prob1NoZero .* log2(prob1NoZero));
            
            % repeat for the second image
            indNoZero2 = histImage2 ~= 0;
            prob2NoZero = histImage2(indNoZero2) / numel(histImage2);
            entropy2 = -sum(prob2NoZero .* log2(prob2NoZero));
            
            % now compute mutual information
            mutualInfo = (entropy1 + entropy2) - jointEntropy;
 
            % assign to output
            out(ii) = mutualInfo;
            
        % angle of intersection (subset the tensor)    
        case 'angle'
            
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
            
            % other measures            
            % compute the angle between connections
            %theta = acos(val);
            % compute cosine distance between angle
            %csdst = pdist2(matm1', matm2', 'cosine'); % cosine distance
            
        case 'tpnorm'
            
            % grab path voxels and values for the first connection
            tp1 = eval([ 'pconn{' num2str(conn1) '}.' label '.profile.' dtype ]);
            tp2 = eval([ 'pconn{' num2str(conn2) '}.' label '.profile.' dtype ]);
            
            % compute the norm between profiles
            tpn = norm(tp1 - tp2);
            
            % other measures
            %tpn = corr(tp1, tp2);               % correlation
            %tpn = mae(tp1 - tp2)                % mean of absolute error
            %tpn = norm(tp1 - tp2, 1)            % sum of absolute error
            %tpn = norm(tp1 - tp2, inf)          % max of absolute error 
            %tpn = dot(tp1, tp2)                 % scalar product
            %tpn = mse(tp1, tp2)                 % mean of square error
            %tpn = sse(tp1 - tp2)                % sum of square error
            %tpn = sqrt(mean((tp1 - tp2) .^ 2)); % RMSE
            %tpn = pdist2(tp1', tp2', 'cosine'); % cosine distance
            
            % combine and fit coefficients to the lines
            % how could these be compared?
            %[ coeff, errs ] = lpc([ tp1, tp2 ], 4);
            
            % store output
            out(ii) = tpn;
            
        otherwise
            
            error('Invalid metric between edges requested.');
            
    end
    
end
time = toc;

disp(['Computed link network metrics in ' num2str(time/60) ' minutes.' ]);

disp('Creating link network graph...');

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

