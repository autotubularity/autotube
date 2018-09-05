function imFrangi = doFrangi(im, opts)
% doFrangi: this function applies the Frangi filter on a given input image.
% 
% inputs:
%       im: a 2D gray-level image.
%       opts: a matlab structure containing the parameters for the frangi
%               filter.
%
% output:
%       imFrangi: frangi-vesselness response.
%
% references: for more details on the method, refer to:
%             Frangi et al. Multiscale Vessel Enhancement Filtering

    if ~isfield(opts, 'tubularity')
        opts.tubularity = [ ];
    end
    
    if ~isfield(opts.tubularity, 'beta1')
        opts.tubularity.beta1 = 0.5;
    end
    
    if ~isfield(opts.tubularity, 'beta2')
        opts.tubularity.beta2 = 15;
    end
    
    if ~isfield(opts.tubularity, 'blackwhite')
        opts.tubularity.blackwhite = false;
    end
    
    if ~isfield(opts.tubularity, 'sigmaInit')
        opts.tubularity.sigmaInit = 1;
    end
    
    if ~isfield(opts.tubularity, 'sigmaRatio')
        opts.tubularity.sigmaRatio = 2; %sqrt(2);
    end
    
    if ~isfield(opts.tubularity, 'sigmaEnd')
        opts.tubularity.sigmaEnd = 5;
    end

    sigmas     = opts.tubularity.sigmaInit:opts.tubularity.sigmaRatio:opts.tubularity.sigmaEnd;
    beta       = 2*(opts.tubularity.beta1)^2;
    correction = 2*(opts.tubularity.beta2)^2;
    
    filtImArr = zeros([size(im) numel(sigmas)]);
    
    for ii=1:numel(sigmas)
        
        [ixx, ixy, iyy] = do_hessian(im, sigmas(ii));
        
        ixx = sigmas(ii)^2*ixx;
        ixy = sigmas(ii)^2*ixy;
        iyy = sigmas(ii)^2*iyy;
        
        [lambda2, lambda1] = eig_image(ixx, ixy, iyy);
        
        % Compute some similarity measures
        lambda1(lambda1==0) = eps;
        Rb = (lambda2./lambda1).^2;
        S2 = lambda1.^2 + lambda2.^2;

        % Compute the output image
        imFilt = exp(-Rb/beta) .*(ones(size(im))-exp(-S2/correction));

        if(opts.tubularity.blackwhite)
            imFilt(lambda1 < 0)=0;
        else
            imFilt(lambda1 > 0)=0;
        end
        
        filtImArr(:,:,ii) = imFilt;
    end
    
    imFrangi = max(filtImArr, [], 3);
end

function [ixx, ixy, iyy] = do_hessian(im, sigma)
% Hessian computation

    hsize = 2*ceil(3*sigma) + 1;
    gf    = fspecial('gauss', hsize, sigma);
    
    [gx, gy]   = gradient(gf);
    [gxx, gxy] = gradient(gx);
    [~, gyy] = gradient(gy);
    
    ixx = imfilter(im, gxx, 'symmetric');
    ixy = imfilter(im, gxy, 'symmetric');
    iyy = imfilter(im, gyy, 'symmetric');
    
end

function [lambda1, lambda2, varargout] = eig_image(ixx, ixy, iyy)
%     | ixx  ixy |
% H = |          |
%     | iyx  iyy |
    
    discriminant = sqrt( (ixx-iyy).^2 + 4*ixy.^2 );
    
    % computing eigenvalues from Hessian Matrix
    lambda1 = 0.5*((ixx + iyy) + discriminant);
    lambda2 = 0.5*((ixx + iyy) - discriminant);
    
    % computing eigenvectors and normalising them.
    v2x = 2*ixy;
    v2y = iyy - ixx + discriminant;
    vmag = sqrt(v2x.^2 + v2y.^2);
    v2x = v2x./(vmag+eps);
    v2y = v2y./(vmag+eps);

    % eigenvectors shall be orthogonal
    eigvect1 = -v2y; 
    eigvect2 = v2x;
    
    % sort eigenvalues and eigenvectors
    check = abs(lambda1)>abs(lambda2);
    lambda1(check) = lambda2(check);
    lambda2(check) = lambda1(check);
    
    eigvect1(check)=v2x(check);
    eigvect2(check)=v2y(check);

    varargout{1} = eigvect1;
    varargout{2} = eigvect2;

end