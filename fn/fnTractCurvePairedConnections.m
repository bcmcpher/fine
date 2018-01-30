function [ pconn, tcrve, pcsf ] = fnTractCurvePairedConnections(fg, pconn, label, nnodes, minNum, clobber)
%feTractShapePairedConnections creates a 
%
% INPUTS:
%     fg      - fiber group in acpc space
%     pconn   - paired connection object created with fg
%     label   - string indicating the fiber groups for which to create tract profiles
%               either:
%                      'all' for all assigned streamlines
%                      'nzw' for non-zero weighted fibers returned by LiFE
%                      'zrw' for all zero weighted fibers removed by LiFE
%     nnodes  - the number of nodes to which tract profiles will be sampled
%     minNum  - the minimum number of streamlines for a profile to be computed.
%               the defualt is 4.
%     clobber - overwrite existing profile label if it already exists
%             
% OUTPUTS:
%     pconn - is the paired connections object with the tract profile(s) added
%             in a field called 'profile'.
%
%     tcrve - debugging output; the cell array that is added internally to pconn
%
%     pcsf  - debugging output: the cell array of all superfiber structures
%
% EXAMPLE:
%
% % load data
% parc  = niftiRead('labels.nii.gz');
% favol = niftiRead('fa.nii.gz');
% fg        = feGet(fe, 'fibers acpc');
% fibers    = fg.fibers;
% fibLength = fefgGet(fg, 'length');
% weights   = feGet(fe, 'fiberweights');
%
% % assign streamlines to edges
% [ pconn, rois ] = feCreatePairedConnections(parc, fibers, fibLength, weights);
%
% % create a profile for every connection with non-zero weighted fibers
% pconn = feTractShapePairedConnections(fg, pconn, 'nzw', 100, 3);
%
% Brent McPherson (c), 2017 - Indiana University
%

% parse optional arguments
if(~exist('nnodes', 'var') || isempty(nnodes))
    nnodes = 100;
end

if(~exist('minNum', 'var') || isempty(minNum))
    minNum = 3;
end

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

if(clobber == 1)
    display(['Overwriting existing shape profiles...']);
end

% check if the first profile label already exists alongside clobber for a
% faster exit from function if things shouldn't be recomputed.
if isfield(pconn{1}.(label), 'profile')
    if isfield(pconn{1}.(label).profile, 'shape') && clobber == 0
        error('Tract shapes for this label already exist. Please set clobber to 1 to explicitly overwrite existing shapes and superfibers.');
    end
end

% check if average length is long enough for curvature values to be reasonable?

% preallocate tract profile output
tcrve = cell(length(pconn), 1);
pcsf = cell(length(pconn), 1);
tccnt = 0;
tctry = 0;

display(['Computing tract curves on ' num2str(length(pconn)) ' connections...']);

tic;    
fibers = fg.fibers;
parfor ii = 1:length(pconn)
    
    % pull field requested
    tmp = pconn{ii}.(label);
    
    % in testing, need minimum of 4 streamlines to compute profile
    if size(tmp.indices, 1) > minNum
                
        % create tract-wise fg
        tfg = fgCreate('fibers', fibers(tmp.indices));
        tfg = dtiReorientFibers(tfg, nnodes);
        
        % create superfiber representation
        pcsf{ii} = dtiComputeSuperFiberRepresentation(tfg, [], nnodes);
        sf_name = [ 'sf_' num2str(pconn{ii}.roi1) '_to_' num2str(pconn{ii}.roi2) ];
        pcsf{ii}.name = sf_name;
        
        % if the variance of the streamlines distance is far, too many
        % streamlines are dropped to reliably compute curves
        
        try
            
            % preallocate outputs
            numfibers = length(tfg.fibers);
            tan  = zeros(nnodes, 3, numfibers);
            norm = zeros(nnodes, 3, numfibers);
            bnrm = zeros(nnodes, 3, numfibers);
            curv = zeros(nnodes, numfibers);
            tors = zeros(nnodes, numfibers);
            weights = zeros(nnodes, numfibers);
            
            % calculate curvature and torsion for each node on each fiber
            % frenet2 embedded at the bottom
            for jj = 1:numfibers
                [ tan(:, :, jj), norm(:,:,jj), bnrm(:,:,jj), ...
                  curv(:, jj),tors(:, jj) ] = ...
                frenet2(tfg.fibers{jj}(1,:)',tfg.fibers{jj}(2,:)',tfg.fibers{jj}(3,:)');
            end
            
            % Each fiber is represented by numberOfNodes, so can easily loop over 1st
            % node, 2nd, etc...
            fc = horzcat(tfg.fibers{:})';
            
            % calculate the weights to be applied to each fibers curvature estimates.
            % these weights are based on the fiber's mahalanobis distance from the
            % fiber core. they will be used to weight that fibers contribution to the
            % core measurement
            for node=1:nnodes
                
                % compute gaussian weights y = mvnpdf(X,mu,SIGMA);
                % returns the density of the multivariate normal distribution with zero
                % mean and identity covariance matrix, evaluated at each row of X.
                X = fc((1:nnodes:numfibers*nnodes)+(node-1), :);
                
                % get the covariance structure from the super fiber group
                fcov = pcsf{ii}.fibervarcovs{1};
                sigma = [fcov(1:3, node)'; ...
                         0 fcov(4:5, node)';...
                         0 0 fcov(6, node)'];
                sigma = sigma + sigma' - diag(diag(sigma));
                mu    = pcsf{ii}.fibers{1}(:, node)';
                
                % weights for the given node.
                weights(node, :) = mvnpdf(X, mu, sigma')';
                
            end
            
            % the weights for each node should add up to one across fibers.
            weightsNormalized = weights./(repmat(sum(weights, 2), [1 numfibers]));
            
            % weight each curvature and torsion calculation by its distance from the core
            tcrve{ii}.curv = sum(curv .* weightsNormalized, 2);
            tcrve{ii}.tors = sum(tors .* weightsNormalized, 2);
            
            % apply weights to curve vectors
            for vec = 1:3 % for x / y / z vector component                
                tan(:,vec,:)  = squeeze(tan(:,vec,:)) .* weightsNormalized;
                norm(:,vec,:) = squeeze(norm(:,vec,:)) .* weightsNormalized;
                bnrm(:,vec,:) = squeeze(bnrm(:,vec,:)) .* weightsNormalized;
            end
            
            % compute the weighted average vector curves
            tcrve{ii}.tan  = mean(tan, 3);
            tcrve{ii}.norm = mean(norm, 3);
            tcrve{ii}.bnrm = mean(bnrm, 3);
            
            % track how many connections are profiled
            tccnt = tccnt + 1;
            tctry = tctry + 1;
            
        catch
            
            warning(['Connection: ' num2str(ii) ' failed to compute curve data.']);
            tctry = tctry + 1;
            tcrve{ii} = nan(nnodes, 1);

        end
        
    else
        
        % skip empty connection
        sf_name = [ 'sf_' num2str(pconn{ii}.roi1) '_to_' num2str(pconn{ii}.roi2) ];
        pcsf{ii} = struct('name', sf_name, 'n', 0, 'fibers', nan(3, nnodes), 'fibervarcovs', nan(6, nnodes));
        tcrve{ii} = nan(nnodes, 1);
        continue
        
    end
    
end
time = toc;

display(['Computed ' num2str(tccnt)  ' of ' num2str(tctry) ' possible tract curves in ' num2str(round(time)/60) ' minutes.']);

clear ii time

%% add tract profile to pconn object

display('Adding tract curves to pconn...');

for ii = 1:length(pconn)
    
    % pull subset field
    tmp = pconn{ii}.(label);
    
    % look for an existing profile field
    if isfield(tmp, 'profile') % if there is one
        
        % pull profile field so what's there isn't lost
        prof = pconn{ii}.(label).profile;
                    
            % clobber is already set and this fxn is the only one to make
            % .shape data field, so don't worry about preserving what's
            % here (for now)
            
            try    
                
                % create the shape field and add the data
                crve = struct('curv', tcrve{ii}.curv, 'tors', tcrve{ii}.tors, ...
                    'tan', tcrve{ii}.tan, 'norm', tcrve{ii}.norm, 'bnrm', tcrve{ii}.bnrm);
                
            catch
                
                % if it's empty, fill in NaNs
                crve = struct('curv', nan(nnodes, 1), 'tors', nan(nnodes, 1), ...
                    'tan', nan(nnodes, 3), 'norm', nan(nnodes, 3), 'bnrm', nan(nnodes, 3));
                
            end
            
        %end
        
        % add tract profile(s) to tmp
        tmp.profile.shape = crve;
        
    % if profile field doesn't exist
    else 
        
        try
            
            % create the shape field and add the data
            crve = struct('curv', tcrve{ii}.curv, 'tors', tcrve{ii}.tors, ...
                'tan', tcrve{ii}.tan, 'norm', tcrve{ii}.norm, 'bnrm', tcrve{ii}.bnrm);
            
        catch
            
            % if it's empty, fill in NaNs
            crve = struct('curv', nan(nnodes, 1), 'tors', nan(nnodes, 1), ...
                'tan', nan(nnodes, 3), 'norm', nan(nnodes, 3), 'bnrm', nan(nnodes, 3));
            
        end
        
        % add profiles and shape with crve sored
        tmp.profile.shape = crve;
        
    end
    
    % reassign superfiber and tmp to a paired connection cell array
    pconn{ii}.(label).suberfiber = pcsf{ii};
    pconn{ii}.(label).profile = tmp.profile;
    
end

end

%% internal functions

function [T,N,B,k,t] = frenet2(x,y,z)
% FRENET - Frenet-Serret Space Curve Invarients
%   
%   [T,N,B,k,t] = frenet2(x,y);
%   [T,N,B,k,t] = frenet2(x,y,z);
% 
%   Returns the 3 vector and 2 scaler invarients of a space curve defined
%   by vectors x,y and z.  If z is omitted then the curve is only a 2D,
%   but the equations are still valid.
% 
%    _    r'
%    T = ----  (Tangent)
%        |r'|
% 
%    _    T'
%    N = ----  (Normal)
%        |T'|
%    _   _   _
%    B = T x N (Binormal)
% 
%    k = |T'|  (Curvature)
% 
%    t = dot(-B',N) (Torsion)
% 
% 
%    Example:
%    theta = 2*pi*linspace(0,2,100);
%    x = cos(theta);
%    y = sin(theta);
%    z = theta/(2*pi);
%    [T,N,B,k,t] = frenet2(x,y,z);
%    line(x,y,z), hold on
%    quiver3(x,y,z,T(:,1),T(:,2),T(:,3),'color','r')
%    quiver3(x,y,z,N(:,1),N(:,2),N(:,3),'color','g')
%    quiver3(x,y,z,B(:,1),B(:,2),B(:,3),'color','b')
%    legend('Curve','Tangent','Normal','Binormal')
% 
% 
% See also: GRADIENT

if nargin == 2
    z = zeros(size(x));
end

% CONVERT TO COLUMN VECTOR
x = x(:);
y = y(:);
z = z(:);

% SPEED OF CURVE
dx = gradient(x);
dy = gradient(y);
dz = gradient(z);
dr = [dx dy dz];

ddx = gradient(dx);
ddy = gradient(dy);
ddz = gradient(dz);
ddr = [ddx ddy ddz];



% TANGENT
T = dr./mag(dr,3);


% DERIVIATIVE OF TANGENT
dTx =  gradient(T(:,1));
dTy =  gradient(T(:,2));
dTz =  gradient(T(:,3));

dT = [dTx dTy dTz];


% NORMAL
N = dT./mag(dT,3);
% BINORMAL
B = cross(T,N);
% CURVATURE
% k = mag(dT,1);
k = mag(cross(dr,ddr),1)./((mag(dr,1)).^3);
% TORSION
%t = dot(-B,N,2);

dddx = gradient(ddx); 
dddy = gradient(ddy); 
dddz = gradient(ddz); 
dddr = [dddx dddy dddz];

t = vdot(cross(dr, ddr), dddr) ./ mag(cross(dr, ddr),1).^2;

end

function N = vdot(A, B)
%row-wise dot-product of A and B
N=zeros(size(A,1),1);
for i=1:size(A,1)
    N(i) = dot(A(i,:), B(i,:));
end
end

function N = mag(T,n)
% MAGNATUDE OF A VECTOR (Nx3)
%  M = mag(U)
N = sum(abs(T).^2,2).^(1/2);
d = find(N==0);
N(d) = eps*ones(size(d));
N = N(:,ones(n,1));
end
