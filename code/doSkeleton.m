function [im_skel, im_branches, varargout] = doSkeleton(im, opts)
% doSkeleton: this function extracts the skeleton from a binary input
% image. Short ramifications can also be pruned.
% 
% inputs:
%       im: a 2D binary image input image containing vessels as forebround objects.
%       opts: a matlab structure containing the parameters for the ramification-prunning:
%           -> opts.skeleton.cleanSpur: 1 to prune ramifications.
%           -> opts.skeleton.spurLength: length of ramifications to be pruned in pixel units.
%
% output:
%       im_skel: 2D binary output image with extracted skeleton.
%       im_branches: binary output image with branches as white pixels.

    if ~isfield(opts, 'skeleton')
        opts.skeleton = [];
    end
    
    if  ~isfield(opts.skeleton, 'cleanSpur')
        opts.skeleton.cleanSpur = 1;
    end

	im_skel = bwmorph(im, 'thin', 'inf');
            
	im_branches = bwmorph(im_skel, 'branchpoints');
    
    if opts.skeleton.cleanSpur == 1
        [im_skel_d, ~]  = getPrunnedBranches(im_skel, opts);
        im_branches_d   = bwmorph(im_skel_d, 'branchpoints');
    else
        im_skel_d = im_skel;
        im_branches_d = im_branches;
    end
    
    varargout{1}  = im_skel_d;
    varargout{2}  = im_branches_d;
end

function [skelD, im_branches] = getPrunnedBranches(im_skel, opts)

    if ~isfield(opts.skeleton, 'spurLength')
        opts.skeleton.spurLength = numel(im_skel);
    end

    B = bwmorph(im_skel, 'branchpoints');
    E = bwmorph(im_skel, 'endpoints');
    
    [rows,cols] = find(E);
    B_loc = (B);

    Dmask = false(size(im_skel));
    for k = 1:numel(cols)
        D = bwdistgeodesic(im_skel,cols(k),rows(k));
        distanceToBranchPt = min(D(B_loc));
        
        if ~opts.skeleton.spurLength 
            Dmask(D < distanceToBranchPt) = true;
        else
            Dmask(D < distanceToBranchPt & distanceToBranchPt < opts.skeleton.spurLength) = true;
        end
    end
    skelD = im_skel - Dmask;
    skelD = bwmorph(skelD,'thin'); 
    
    [yy,xx] = find(B); 
    im_branches = zeros(size(im_skel));
    im_branches(sub2ind(size(im_branches), yy, xx)) = 1;
    im_branches = logical(im_branches);

end