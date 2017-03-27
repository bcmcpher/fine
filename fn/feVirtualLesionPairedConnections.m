function [ pconn, vlout ] = feVirtualLesionPairedConnections(fe, pconn, label)
%feVirtualLesionPairedConnections runs virtual lesion on a field of pconn indices.
%
% INPUTS:
% - 'fe' is a fit fe structure that is used to build pconn
% - 'pconn' is the paired connection structure that stores fiber indices
%   between each unique pair of connections in the 
% - 'label' is the set of indices in pconn to perform virtual lesions for
%
% OUTPUTS:
% - 'pconn' is the paired connections object with VL data stored and added
%   to the matrix field for generating networks
% - 'vlout' is the virtual lesion output only as a cell array. For debugging.
%

% run parallelized virtural lesion

% preallocate virtual lesion output
vlout = cell(length(pconn), 1);
vlcnt = 0;

tic;
parfor ii = 1:length(pconn)
    
    % pull field requested
    tmp = getfield(pconn{ii}, label);
    
    if sum(tmp.weights) == 0
        
        % set to zero and continue
        %vlout{ii} = {};
        vlout{ii}.s.mean  = 0;
        vlout{ii}.em.mean = 0;
        vlout{ii}.j.mean  = 0;
        vlout{ii}.kl.mean = 0;
        continue
        
    else
        % CHECK ORDER OF weVL / ewoVL WITH CESAR
        % compute a virtual lesion
        [ ewVL, ewoVL ] = feComputeVirtualLesion(fe, tmp.indices);
        vlout{ii} = feComputeEvidence(ewoVL, ewVL);
        vlcnt = vlcnt + 1;
    end
    
end
time = toc;

display(['Computed ' num2str(vlcnt) ' vitual lesions in ' num2str(round(time)/60) ' minutes.']);

clear ii vlcnt ewVL ewoVL time

% add virtual lesion to pconn object
for ii = 1:length(pconn)
    
    % pull subset field
    tmp = getfield(pconn{ii}, label);
    
    % add whole vl output
    tmp.vl = vlout{ii};
    
    % assign virtual lesion matrix entries
    tmp.matrix.soe = vlout{ii}.s.mean;
    tmp.matrix.emd = vlout{ii}.em.mean;
    tmp.matrix.kl = mean(vlout{ii}.kl.mean);
    tmp.matrix.jd = mean(vlout{ii}.j.mean);
    
    % reassign virtual lesions to paired connection array
    pconn{ii} = setfield(pconn{ii}, label, tmp);
    
end

end
