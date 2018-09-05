function im_den = doDenoising(im, opts)
% doDenoising: this function denoises an input image, i.e. reduces noise.
% By default a BM3D-based denoising is applied.
% 
% inputs:
%       im: a 2D gray-level image.
%       opts: a matlab structure containing the type of denoising to be
%               applied.
%       -> opts.denoising.type: (string) type of denoising to be applied: 'BM3D' or 'wiener' filter.
%
% output:
%       im_den: denoised image.

    if ~isfield(opts, 'denoising')
       opts.denoising= []; 
    end

    if ~isfield(opts.denoising, 'type')
       opts.denoising.type = 'BM3D'; 
    end
    
    if ~isfield(opts.denoising, 'bm3d')
       opts.denoising.bm3d = ''; 
    end
    
    if ~isfield(opts.denoising.bm3d, 'sigma')
       opts.denoising.bm3d.sigma = 100; 
    end
    
    if strcmpi(opts.denoising.type, 'BM3D')
        im = im2double(im);
        if ispc
            [~, im_den] = BM3D_win(1, im, opts.denoising.bm3d.sigma);
        elseif ismac || islinux
            [im_den, ~] = BM3D_mac_linux(im, opts.denoising.bm3d.sigma/255.0, 'np');
        else
            error('Operating System not supported!');
        end
    elseif strcmpi(opts.denoising.type, 'wiener')
        [im_den,~] = wiener2(im, [5 5]);
    else
        error('Unkown denoising method: %s\n', opts.type);
    end

end