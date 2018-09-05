function im_illu = doIlluCorrection(im, opts)
% doIlluCorrection: this function normalizes illumination from a gray-level image
% using the top-hat transform.
% 
% inputs:
%       im: a 2D gray-level image.
%       opts: a matlab structure containing the parameters for the
%             illumination correction. The default size of the circular
%             kernel is of 51.
% output:
%       im_illu: 2D illuminated-normalized image.

    if ~isfield(opts, 'illumination')
        opts.illumination = [];
    end

    if ~isfield(opts.illumination,'seSzeIllu')
        se = strel('disk', 51);
    else
        se = strel('disk', opts.illumination.seSzeIllu);
    end
    
	op  = imopen(im, se);
    im_illu = im - op;
    
end