function im_hull = doHull(im, opts)
% doHull: this function computes the convex hull from a binary image.
% 
% inputs:
%       im: a 2D gray-level image.
%       opts: a matlab structure containing the parameters for convex-hull
%             method. It can be localy computed or globally. By default the
%             global convex-hull is used.
%           -> opts.hull.type (string) convex-hull type 'local'/'global'.
%
% output:
%       im_hull: convex-hull of binary image.


    if ~isfield(opts, 'hull')
        opts.hull = [];
    end
    
    if ~isfield(opts.hull, 'type')
        opts.hull.type = 'global';
    end
    
    if strcmpi(opts.hull.type, 'global')
        im_hull = bwconvhull(im);
    elseif strcmpi(opts.hull.type, 'local')
        im_hull = bwconvhull(im, 'objects');
    else
        error('Unknown convex-hull type: %s!\n', opts.hull.type);
    end
    
end