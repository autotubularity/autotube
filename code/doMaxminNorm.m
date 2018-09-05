function im_out = doMaxminNorm(im_in, newmin, newmax, varargin)
% doMaxminNorm: this function normalizes a given 2D image into the range [newmin, newmax]
% 
% inputs:
%       im: a 2D gray-level image.
%       newmin: new minimun pixel value.
%       newmax: new maximum pixel value.
%
% output:
%       im_out: 2D max-min normalized image.

    im_in = single(im_in);
    numvarargs = length(varargin);
    if numvarargs > 2
        error('do_maxmin_norm:TooManyInputs', ...
        'requires at most 2 optional inputs');
    end
    
    % set defaults for optional inputs
    optargs = {0 255};
    optargs(1:numvarargs) = varargin;
    [defmin, defmax] = optargs{:};

    if numvarargs
        minval = defmin;
        maxval = defmax;  
    else
        minval = min(im_in(:));
        maxval = max(im_in(:));      
    end
    
    im_out = (im_in-minval)*((newmax-newmin)/(maxval-minval)) + newmin;
    
end