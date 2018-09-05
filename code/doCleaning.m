function im_clean = doCleaning(im, opts)
% doCleaning: this function post-process a binary image by removing small
%       objects or by filling holes.
% 
% inputs:
%       im: a 2D binary input image.
%       opts.cleaning: a matlab structure containing some parameters for the post-processing.
%       -> opts.cleaning.minAreaPercent: minimal percenta size of the
%                   total image size for an object to not be removed.
%       -> opts.cleaning.minHoleSize: minimum hole size in pixels to be
%                   filled.
%
% output:
%       im_clean: post-processed binary image.

    if ~isfield(opts, 'cleaning')
        opts.cleaning = [];
    end

    if ~isfield(opts.cleaning, 'minAreaPercent')
        opts.cleaning.minAreaPercent = 0.01/100;
    end
    
    if ~isfield(opts.cleaning, 'minHoleSize')
        opts.cleaning.minHoleSize = 0;
    end

    bw          = logical(im);
    stats       = regionprops(bw,'area','Perimeter');
    
    cc          = bwconncomp(bw, 8);
    L           = labelmatrix(cc);
    
    allAreas    = [stats.Area];
    
    minArea     = opts.cleaning.minAreaPercent * numel(im);

    keepMask    = allAreas > minArea;
    
    bwSmallSegs = ismember(L, find(keepMask));
    
    leftIdx     = unique(L(:,1));   leftIdx(leftIdx==0) = [];
    rightIdx    = unique(L(:,end)); rightIdx(rightIdx==0) = [];
    topIdx      = unique(L(1,:));   topIdx(topIdx==0) = [];
    bottomIdx   = unique(L(end,:)); bottomIdx(bottomIdx==0) = [];
    
    bwSmallSegs(ismember(L, leftIdx)) = 1;
    bwSmallSegs(ismember(L, rightIdx)) = 1;
    bwSmallSegs(ismember(L, topIdx)) = 1;
    bwSmallSegs(ismember(L, bottomIdx)) = 1;
    
    im_clean      = (bw & bwSmallSegs);
    
    if opts.cleaning.minHoleSize 
        minimum_hole_size   = opts.cleaning.minHoleSize;
        allfilled           = imfill(im_clean, 'holes');
        allholes            = allfilled & ~im_clean;
        threshold_holes     = bwareaopen(allholes, minimum_hole_size);
        small_holes         = allholes & ~threshold_holes;
        im_clean              = im_clean | small_holes;
    end
    
end