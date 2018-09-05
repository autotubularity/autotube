function im_adj = doAdjust(im, opts)
% doAdjust: this function adjust the intensity of a given input image. By
% default an "adaptive histogram equalization" transform is applied.
% 
% inputs:
%       im: a 2D gray-level image.
%       opts: a matlab structure containing the type of intensity
%           -> opts.adjustment (string) adjustment to be applied, 'autocontrast', 'global', or 'adaptive' (default)
%
% output:
%       im_adj: intensity-adjusted image.

    if ~isfield(opts, 'adjustment')
       opts.adjustment = 'autocontrast'; 
    end
    
    opts.adjustment = lower(opts.adjustment);
    
    if isempty(strfind(opts.adjustment, 'autocontrast')) == 0
        im_adj = imadjust(im); 
    elseif isempty(strfind(opts.adjustment, 'global')) == 0
        im_adj = histeq(im);
    elseif isempty(strfind(opts.adjustment, 'adaptive')) == 0
        im_adj = adapthisteq(im);
    else
        im_adj = adapthisteq(im);
    end
  
end