function [im_bw, varargout] = doThresh(im, opts)
% doThresh: this function performs image thresholding on a gray-level
% image.
% 
% inputs:
%       im: a 2D gray-level image.
%       opts: a matlab structure containing parameters for the thresholding
%             type:
%          -> opts.thresh: otsu/multi-otsu/kittler/adaptive. By default otsu is used.
%
% output:
%       im_bw: 2D binarized/thresholded image.

    if ~isfield(opts, 'thresh')
       opts.thresh = 'kittler'; 
    end
    
    if ~isfield(opts, 'fillSmallHoles')
       opts.fillSmallHoles = true; 
    end
    
    if ~isa(im,'uint8')
        im = im2uint8(im);
    end
    
    opts.thresh = lower(opts.thresh);

    im_gray = [];
    if strcmpi(opts.thresh, 'otsu')
        level = graythresh(im);
        im_bw = imbinarize(im,level);
    elseif strcmpi(opts.thresh, 'multiotsu')
        thresh = multithresh(im, 3);
        seg_I = imquantize(im, thresh);
        im_bw = seg_I > 1;
    elseif strcmpi(opts.thresh, 'kittler')
        im = im2uint8(im);
        [optThresh, ~ ] = kittler( im );
        im_bw       = zeros(size(im));
        im_bw(im > optThresh) = 1;
        im_gray     = im;
        im_gray(im < optThresh) = 0;
    elseif strcmpi(opts.thresh, 'adaptive')
        im = im2uint8(im);
        T = adaptthresh(im, 0.5);
        im_bw = imbinarize(im,T);
    else
        level = graythresh(im);
        im_bw = imbinarize(im,level);
    end
    
    if opts.fillSmallHoles
       filled     = imfill(im_bw, 'holes');
       holes      = filled & ~im_bw;
       bigholes   = bwareaopen(holes, 30);
       smallholes = ~bigholes & holes;
       im_bw      = im_bw | smallholes;
    end
    
    %# smoothing thresholded boundaries.
    im_bw = imopen(im_bw, ones(2,2)); 
    
    varargout{1} = im_gray;
end