function comb = doOverlay(im, mask, alpha, overlaycolor)
% doOverlay: this function does an alpha blending on a gray-level image
% on the regions defined by a binary mask.
% 
% inputs:
%       im: a gray-level image.
%       mask: a binary image containing the regions to be blended.
%       alpha: alpha blending percentage, should be between [0,1].
%       overlaycolor: 'yellow', 'magenta', 'cyan', 'green', 'blue', 'red'
%
% output:
%       comb: alpha-blended image


    DEFAULT_COLOR = [0 1 0];
    if nargin < 4
        overlaycolor = DEFAULT_COLOR;
    elseif ischar(overlaycolor)
        switch overlaycolor
            case {'y','yellow'}
                overlaycolor = [1 1 0];
            case {'m','magenta'}
                overlaycolor = [1 0 1];
            case {'c','cyan'}
                overlaycolor = [0 1 1];
            case {'r','red'}
                overlaycolor = [1 0 0];
            case {'g','green'}
                overlaycolor = [0 1 0];
            case {'b','blue'}
                overlaycolor = [0 0 1];
            case {'w','white'}
                overlaycolor = [1 1 1];
            case {'gr','gray'}
                overlaycolor = [0.8 0.8 0.8];
            case {'k','black'}
                overlaycolor = [0 0 0];
            otherwise
                disp('Unrecognized color specifier; using default.');
                overlaycolor = DEFAULT_COLOR;
        end
    end

    im           = im2double(im);
    mask         = logical(mask);
    overlaycolor = im2double(overlaycolor);
    
    if size(im,3) == 1
       im = cat(3,im,im,im); 
    end
    
	if size(alpha,2) == 1
        alpha = [1.0, alpha];
	end
    
    alphaFr     = alpha(1);
    alphaBkr    = alpha(2);
    
    r_colour    = overlaycolor(1);
    g_colour    = overlaycolor(2);
    b_colour    = overlaycolor(3);
    
    r_channel   = r_colour*mask;
    g_channel   = g_colour*mask;
    b_channel   = b_colour*mask;
    
    mask_3      = cat(3, r_channel, g_channel, b_channel);
    
    comb        = zeros(size(mask_3));
    comb(:,:,1) = im(:,:,1)*alphaFr + mask_3(:,:,1)*(1.0-alphaBkr);
    comb(:,:,2) = im(:,:,2)*alphaFr + mask_3(:,:,2)*(1.0-alphaBkr);
    comb(:,:,3) = im(:,:,3)*alphaFr + mask_3(:,:,3)*(1.0-alphaBkr);
    
end