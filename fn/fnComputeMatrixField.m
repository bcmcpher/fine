function [ out ] = fnComputeMatrixField(conn, label)
%fnComputeMatrixField creates the valuest that are assigned to adjacency
% matrices for the pconn list object. New measures should be added here.
%
% If pvoxels exists, I can compute volume, ADF, any micro-structure average
% How should this fxn do that? 
% 

% calculate combined size of ROI
psz = conn.roi1sz + conn.roi2sz;

% pull the streamline count
cnt = size(conn.(label).indices, 1);

% pull the average streamline length
len = mean(conn.(label).lengths);

% create 1 / sum of lengths for Hagmann's correction
if isempty(conn{ii}.(label).lengths) 
    dln = 0;
else
    dln = sum(1 / conn{ii}.(label).lengths);
end

% compute tract volume - number of voxels
%size(conn.pvoxel, 1)
% scale by voxel dimesions too?

% compute average microstructural value
% needs whole fe structure?
% have to do separately?

% create all streamline counts
out.matrix.count = cnt;
out.matrix.density = (2 * cnt) / psz;
out.matrix.length = len;
out.matrix.denlen = (2 / psz) * dln;
% out.matrix.afdconn = []; % tract volume / len
% out.matrix.avgmsv = []; % averaged microstructural value

end

