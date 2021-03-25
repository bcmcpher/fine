function [ netw ] = fnComputeMatrixEdges(netw, den, life, clobber)
%fnComputeMatrixField creates the values that are assigned to adjacency
% matrices for the pconn list object. New measures should be added here.
%

%% parse arguments

if(~exist('clobber', 'var') || isempty(clobber))
    clobber = 0;
end

if(~exist('den', 'var') || isempty(den))
    den = 'vol';
end

if(~exist('life', 'var') || isempty(life))
    life = [];
end

if clobber == 0 && isfield(netw.pconn{1}, 'matrix')
    error('Matrix field has already been computed. Set clobber to true to recompute values.');
end

if clobber == 1 && isfield(netw.pconn{1}, 'matrix')
    disp('Matrix field will be overwritten for every edge...');
end

%% compute the matrix fields

% grab the sorted indice pairs
pairs = netw.parc.pairs;
nconn = size(netw.edges, 1);

% for every connection
for conn = 1:nconn
    
    % pull the rois for the pair and the edge information
    roi1 = netw.nodes{pairs(conn, 1)};
    roi2 = netw.nodes{pairs(conn, 2)};
    edge = netw.edges{conn};
    
    % calculate combined size of ROI (voxel count / volume in mm^2)
    psz = roi1.size + roi2.size;
    pvl = roi1.volume + roi2.volume;
    
    % pull the streamline count
    cnt = size(edge.fibers.indices, 1);
    
    % pull the average streamline length
    len = mean(edge.fibers.lengths, 'omitnan');
    
    % compute standard error of length
    len_se = std(edge.fibers.lengths, 'omitnan') / sqrt(size(edge.fibers.lengths, 1));
    
    % Replace mean/se of length w/ 0 if it's empty
    if isnan(len)
        len = 0;
        len_se = 0;
    end
        
    % create 1 / sum of lengths for Hagmann's correction
    if isempty(edge.fibers.lengths)
        dln = 0;
    else
        dln = sum(1 / edge.fibers.lengths);
    end

    % create all count edge measures
    matrix.count = cnt;
    matrix.length = len;
    matrix.len_se = len_se;

    % streamline density normalized by volume by default, optionally voxel count
    switch den
        case {'sz', 'size'}
            % normalize streamline count by number of voxels in both nodes
            matrix.density = (2 * cnt) / psz;
            matrix.denlen = (2 / psz) * dln;
        otherwise
            % normalize streamline count by computed volume of both nodes
            matrix.density = (2 * cnt) / pvl;
            matrix.denlen = (2 / pvl) * dln;
    end
    
    %% deal with LiFE weights being defined and passed
    
    % if a field is passed for LiFE weights, estimate the weight based measures
    if ~isempty(life)
        
        % logically index for LiFE weights above zero
        wght = edge.fibers.(life) > 0;
                
        % pull the streamline count
        lcnt = size(edge.fibers.indices(wght), 1);
    
        % pull the average streamline length
        llen = mean(edge.fibers.lengths(wght), 'omitnan');
        
        % compute standard error of length
        llen_se = std(edge.fibers.lengths(wght), 'omitnan') / sqrt(size(edge.fibers.lengths(wght), 1));
        
        % Replace mean/se of length w/ 0 if it's empty
        if isnan(len)
            llen = 0;
            llen_se = 0;
        end
        
        % create 1 / sum of lengths for Hagmann's correction
        if isempty(edge.fibers.lengths(wght))
            ldln = 0;
        else
            ldln = sum(1 / edge.fibers.lengths(wght));
        end
        
        % create all count edge measures
        matrix.life_count = lcnt;
        matrix.life_length = llen;
        matrix.life_len_se = llen_se;
        
        % streamline density normalized by volume by default, optionally voxel count
        switch den
            case {'sz', 'size'}
                % normalize streamline count by number of voxels in both nodes
                matrix.life_density = (2 * lcnt) / psz;
                matrix.life_denlen = (2 / psz) * ldln;
            otherwise
                % normalize streamline count by computed volume of both nodes
                matrix.life_density = (2 * lcnt) / pvl;
                matrix.life_denlen = (2 / pvl) * ldln;
        end
        
        % unique to LiFE weight measures used in Qi et al. (2016)
        matrix.life_fc = sum(edge.fibers.(life), 'omitnan');
        matrix.life_fw = mean(edge.fibers.(life), 'omitnan');
        %matrix.life_fw = std(edge.fibers.(life), 'omitnan') / sqrt(size(edge.fibers.(life), 1));
        
    end
    
    % option to look at removed (zero-weighted) fibers only?
    
    %% deal with virtual lesion being stored
    
    % if norm virtual lesion is present, store the values
    if isfield(edge, 'vl')
        matrix.soe_mn = edge.vl.s.mean;
        matrix.soe_sd = edge.vl.s.std;
        matrix.emd = edge.vl.em.mean;
        matrix.kld = edge.vl.kl.mean;
        matrix.jef = edge.vl.j.mean;
    end
        
    %% if volume is present, compute the volume measures
    
    if isfield(edge, 'volume')
        matrix.volume = edge.volume.volume;
        matrix.prpvol = edge.volume.prpvol;
    end
    
    %% if microstructure is present, compute mn/std microstructural values
    
    % if microstructure central tendencies are stored
    if isfield(edge, 'micros')
        
        % pull all microstructure labels
        vals = fieldnames(edge.micros);
        
        % for every labels
        for val = 1:length(vals)
            
            lab = vals{val};
            
            % store the values in the matrix field
            matrix.([ lab '_mn' ]) = edge.micros.(lab).mn;
            matrix.([ lab '_md' ]) = edge.micros.(lab).md;
            matrix.([ lab '_sd' ]) = edge.micros.(lab).sd;
            
        end
        
    end
    
    %% what else is optional?
    
    % assign the output back
    edge.matrix = matrix;
    
    % reassign edge
    netw.edges{conn} = edge;
    
end

end

