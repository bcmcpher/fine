function [ netw ] = fnTractCurveEdges(netw, fg, nnodes, minNum, clobber)
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
%               the defualt is 3.
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
    disp('Overwriting any existing shape profiles...');
end

% extract connection list and fibers
pconn = netw.edges;
fibers = fg.fibers;

% check if the first profile label already exists alongside clobber for a
% faster exit from function if things shouldn't be recomputed.
if isfield(pconn{1}, 'curve') && clobber == 0
    error('Tract shapes already exist. Please set clobber to 1 to explicitly overwrite existing shapes and superfibers.');
end

% check if profile labels curv and tors exist, error w/o clobber
if isfield(pconn{1}, 'profile') && clobber == 0
    if isfield(pconn{1}.profile, 'curv') || isfield(pconn{1}.profile, 'tors')
        error('Profiles with the label ''curv'' or ''tors'' already exist. Set clobber = 1 to overwrite pre-existing profile(s).');
    end
    disp('Some profiles have already been computed.');
end

% if a labeled field exists and clobber = 1, print a notification
if (isfield(pconn{1}.profile, 'curv') || isfield(pconn{1}.profile, 'tors')) && clobber == 1
    disp('Profiles with the label ''curv'' or ''tors'' are being recomputed.');
end

% check if average length is long enough for curvature values to be reasonable?

% preallocate try counts
tccnt = 0;
tctry = 0;

% pull a quick count of edges w/ streamlines greater than the minimum
pc = sum(cellfun(@(x) size(x.fibers.indices, 1) > minNum, netw.edges));
disp(['Computing tract curves on ' num2str(pc) ' connections...']);

tic;    
for ii = 1:size(pconn, 1)
    
    % pull edge requested
    edge = pconn{ii};
    r1_idx = netw.parc.pairs(ii, 1);
    r2_idx = netw.parc.pairs(ii, 2);
    sf_name = strcat(netw.nodes{r1_idx}.name, '-', netw.nodes{r2_idx}.name);
    
    % as long as enough streamlines are in the edge
    if size(edge.fibers.indices, 1) > minNum
                
        % create edge fg
        tfg = fgCreate('fibers', fibers(edge.fibers.indices));
        [ tfg, epi ] = dtiReorientFibers(tfg, nnodes); % in theory don't have to resample
      
        % pull roi centers in acpc space
        roi1 = netw.nodes{r1_idx}.centroid.acpc';
        roi2 = netw.nodes{r2_idx}.centroid.acpc';
        
        % find the distance between the start of the profile and each roi center
        epi_roi1 = norm(epi - roi1);
        epi_roi2 = norm(epi - roi2);
        
        % if roi2 is closer to the start of the tract profile than roi1, lrflip the fibers
        if (epi_roi2 < epi_roi1)
            tfg.fibers = cellfun(@(x) fliplr(x), tfg.fibers, 'UniformOutput', false);
        end
        
        % create superfiber representation
        if ~isfield(edge, 'superfiber')
            edge.superfiber = dtiComputeSuperFiberRepresentation(tfg, [], nnodes);
            edge.superfiber.name = sf_name;
        end
        
        % if the variance of the streamlines distance is far, too many
        % streamlines are dropped to reliably compute curves
        
        try
            
            % preallocate outputs
            numfibers = size(tfg.fibers, 1);
            tan  = zeros(nnodes, 3, numfibers);
            nrml = zeros(nnodes, 3, numfibers);
            bnrm = zeros(nnodes, 3, numfibers);
            curv = zeros(nnodes, numfibers);
            tors = zeros(nnodes, numfibers);
            weights = zeros(nnodes, numfibers);
            
            % calculate curvature and torsion for each node on each fiber
            % frenet2 embedded at the bottom
            for jj = 1:numfibers
                [ tan(:, :, jj), nrml(:,:,jj), bnrm(:,:,jj), ...
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
                fcov = edge.superfiber.fibervarcovs{1};
                sigma = [fcov(1:3, node)'; ...
                         0 fcov(4:5, node)';...
                         0 0 fcov(6, node)'];
                sigma = sigma + sigma' - diag(diag(sigma));
                mu    = edge.superfiber.fibers{1}(:, node)';
                
                % weights for the given node.
                weights(node, :) = mvnpdf(X, mu, sigma')';
                
            end
            
            % the weights for each node should add up to one across fibers.
            weightsNormalized = weights./(repmat(sum(weights, 2), [1 numfibers]));
            
            % weight each curvature and torsion calculation by its distance from the core
            edge.profile.curv = sum(curv .* weightsNormalized, 2);
            edge.profile.tors = sum(tors .* weightsNormalized, 2);
            
            % apply weights to curve vectors
            for vec = 1:3 % for x / y / z vector component                
                tan(:,vec,:)  = squeeze(tan(:,vec,:)) .* weightsNormalized;
                nrml(:,vec,:) = squeeze(nrml(:,vec,:)) .* weightsNormalized;
                bnrm(:,vec,:) = squeeze(bnrm(:,vec,:)) .* weightsNormalized;
            end
            
            % compute the weighted average vector curves
            edge.curve.tan  = mean(tan, 3);
            edge.curve.norm = mean(nrml, 3);
            edge.curve.bnrm = mean(bnrm, 3);
            
            % track how many connections are profiled
            tccnt = tccnt + 1;
            tctry = tctry + 1;
            
        catch
            
            warning(['Connection: ' num2str(ii) ' failed to compute curve data.']);
            tctry = tctry + 1;
            edge.profile.curv = nan(nnodes, 1);
            edge.profile.tors = nan(nnodes, 1);
            edge.curve.tan = nan(nnodes, 3);
            edge.curve.norm = nan(nnodes, 3)';
            edge.curve.bnrm = nan(nnodes, 3)';
            
        end
        
    else
        
        % too few streamlines exist to compute a profile / it's empty
        if ~isempty(edge.fibers.indices)
            warning([ 'Edge ' num2str(ii) ' is not empty: ' num2str(size(edge.fibers.indices, 1)) ' streamline(s) present; less than ' num2str(minNum) ]);
        end
        
        % skip empty connection, fill in empty values
        %edge.superfiber = struct('name', sf_name, 'n', 0, 'fibers', nan(3, nnodes), 'fibervarcovs', nan(6, nnodes));
        edge.superfiber = [];
        edge.profile.curv = []; %nan(nnodes, 1);
        edge.profile.tors = []; %nan(nnodes, 1);
        edge.curve.tan = []; %nan(nnodes, 3);
        edge.curve.norm = []; %nan(nnodes, 3);
        edge.curve.bnrm = []; %nan(nnodes, 3);
                
    end
    
    % reassign edge with curvature statistics added
    pconn{ii} = edge;
    
end
time = toc;

netw.edges = pconn;

disp(['Computed ' num2str(tccnt)  ' of ' num2str(tctry) ' possible edge curves in ' num2str(round(time)/60) ' minutes.']);

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
