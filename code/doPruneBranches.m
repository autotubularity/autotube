function im_branch_clean = doPruneBranches(im_branches, opts)
% doPruneBranches: this function prunes branch-points from a binary input
% image containing branch-points as white pixels. More precisely,
% branch-points falling within a circular radius are averaged into one single
% new branch-point.
% 
% inputs:
%       im_branches: a 2D binary image containing branch-points as white pixels.
%       opts: a matlab structure containing the parameters for the branch-prunning:
%           -> opts.branches.pruneDist: circular radius for branch-point search in pixel units.
%
% output:
%       im_branch_clean: 2D pruned-branch point binary image.

    if ~isfield(opts, 'branches')
        opts.branches = [];
    end
    
    if ~isfield(opts.branches, 'pruneDist')
        opts.branches.pruneDist = 15;
    end

    if ~islogical(im_branches)
        im_branches = logical(im_branches);
    end
    
    im_branch_clean = zeros(size(im_branches));

    [yy,xx]=find(im_branches==1);
    points = [yy, xx];
    no_points = size(points,1);
    
    [idxs,dist]   = knnsearch(points, points,'k', no_points, 'distance','chebychev'); %'chebychev' 'euclidean'
    dist_bin_near = dist < opts.branches.pruneDist;
    idx_val_near  = (idxs.*dist_bin_near);
    
    triang_mat  = zeros(no_points, no_points);
    visited_pts = ones(1, no_points);
    
    for ii=1:no_points
        neighs      = idx_val_near(ii,:);
        neighs(neighs==0) = [];
        neighs      = neighs(2:end);
        neighs      = neighs(visited_pts(neighs)~=0);
        noNeighs    = numel(neighs);
        keep        = ones(1, noNeighs);
        
        if noNeighs >= 1
           for jj=1:noNeighs
               [knneigh, knneigh_rank]      = find(idx_val_near==neighs(jj));
               [diffneigh, idx_diffneigh]   = setdiff(knneigh, [neighs, ii]);
               if ~isempty(diffneigh)
                   for kk=1:numel(diffneigh)
                        neigh_idx = knneigh_rank(idx_diffneigh(kk));
                        if dist(diffneigh(kk),neigh_idx) < dist(ii, jj+1)
                            keep(jj) = 0;
                        end
                   end
               end
           end
        end
        neighs = neighs(keep==1);
        
        if visited_pts(ii) == 1
            triang_mat(ii,ii) = 1;
        end
        visited_pts(ii) = 0;
        
        if ~isempty(neighs)
            if isrow(neighs), neighs=neighs'; end
            idx = [repmat(ii,numel(neighs), 1), neighs];
            triang_mat(sub2ind(size(triang_mat), idx(:,1), idx(:,2))) = 1;
            visited_pts(neighs) = 0;
        end
        
    end
    
    for ii=1:no_points
        yy_valid = yy(triang_mat(ii,:)~=0);
        xx_valid = xx(triang_mat(ii,:)~=0);
        
        if ~isempty(yy_valid) && ~isempty(xx_valid)
            yy_mean = mean(yy_valid);
            xx_mean = mean(xx_valid);
            im_branch_clean(round(yy_mean), round(xx_mean)) = 1;
        end
    end
    
end