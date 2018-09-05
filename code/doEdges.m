function im_edges = doEdges(im, opts)
% doEdges: this function detects edges on an image as probability maps,
%           i.e. the confidence of an edge being an edge.
% 
% inputs:
%       im: a 2D gray-level image.
%       opts: a matlab structure containing the parameters for the edge detection
%
% output:
%       im_edges: edge probability map.
%
% references: for more details on the method, refer to:
%             Piotr Dollar and C. Lawrence Zitnick: Structured Forests for Fast Edge Detection

    im = (im2uint8(im));

    if isfield(opts, 'model')
        model = opts.model;
    end
    
    % set detection parameters (can set after training)
    model.opts.multiscale=0;          % for top accuracy set multiscale=1
    model.opts.sharpen=2;             % for top speed set sharpen=0
    model.opts.nTreesEval=4;          % for top speed set nTreesEval=1
    model.opts.nThreads=4;            % max number threads for evaluation
    model.opts.nms=0;                 % set to true to enable nms

    if size(im,3)==1
       im = cat(3,im,im,im); 
    end

    tic, im_edges = edgesDetect(im, model); toc

end