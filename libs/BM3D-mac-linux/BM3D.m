function [y_est, y_hat] = BM3D(z, PSD, profile)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  BM3D is an algorithm for attenuation of additive spatially correlated
%  stationary (aka colored) Gaussian noise in grayscale and multichannel images.
%
%
%  FUNCTION INTERFACE:
%
%  y_est = BM3D(z, sigmaPSD, profile)
%
%          'z'  noise image
%
%  INPUT ARGUMENTS:
%
%         'z' : noisy image (M x N double array, intensities in range [0,1])
%  'sigmaPSD' : noise power spectral density (M x N double nonnegative array)
%   'profile' : 'np' --> Normal Profile
%               'lc' --> Fast Profile (slightly lower quality)
%
%  OUTPUT:
%      'y_est'  denoised image  (M x N double array)
%
%
%
%  BASIC SIMULATION EXAMPLES:
%
%     Case 1)
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive colored random noise generated as convution of AWGN against with kernel 'k'
%
%      k=[-1;2;-1]*[1 4 1]/100;   % e.g., a diagonal kernel
%      z=y+imfilter(randn(size(y)),k(end:-1:1,end:-1:1),'circular');
%
%      % define 'sigmaPSD' from the kernel 'k'
%
%      sigmaPSD=abs(fft2(k,size(z,1),size(z,2))).^2*numel(z);
%
%      % Denoise 'z'
%      y_est = BM3D(z, sigmaPSD);
%
%
%     Case 2)
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive colored random noise generated as convution of AWGN against with kernel 'k'
%      [x2, x1]=meshgrid(ceil(-size(y,2)/2):ceil(size(y,2)/2)-1,ceil(-size(y,1)/2):ceil(size(y,1)/2)-1)
%      sigmaPSD=ifftshift(exp(-((x1/size(y,1)).^2+(x2/size(y,2)).^2)*10))*numel(y)/100;
%      z=y+real(ifft2(fft2(randn(size(y))).*sqrt(sigmaPSD)/sqrt(numel(y))));
%
%      % Denoise 'z'
%      y_est = BM3D(z, sigmaPSD);
%
%     Case 3) If 'sigmaPSD' is a singleton, this value is taken as sigma and it is assumed that the noise is white variance sigma^2.
%
%      % Read a grayscale noise-free image
%
%      y=im2double(imread('cameraman.tif'));
%
%      % Generate noisy observations corrupted by additive white Gaussian noise with variance sigma^2
%      sigma=0.1;
%      z=y+sigma*randn(size(y));
%
%      y_est = BM3D(z, sigma);
%
%      % or, equivalently,
%      sigmaPSD = ones(size(z))*sigma^2*numel(z)
%      y_est = BM3D(z, sigmaPSD)
%
%
%      Case 4)   MULTICHANNEL PROCESSING
%
%      y_est = BM3D(cat(3, z1, z2, z3), PSD, 'np'); 
%
%      Multiple PSDs are optionally handled in the same way.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright (c) 2006-2018 Tampere University of Technology.
% All rights reserved.
% This work (software, material, and documentation) shall only
% be used for nonprofit noncommercial purposes.
% Any unauthorized use of this work for commercial or for-profit purposes
% is prohibited.
%
% AUTHORS:
%     M. Mäkinen, L. Azzari, K. Dabov, A. Foi
%     email: alessandro.foi@tut.fi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Quality/complexity trade-off profile selection
%%%%
%%%%  'np' --> Normal Profile (balanced quality)
%%%%  'lc' --> Low Complexity Profile (fast, lower quality)
%%%%

if ~exist('profile','var')
    profile         = 'np'; %% default profile
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Following are the parameters for the Normal Profile.
%%%%

%%%% Select transforms ('dct', 'dst', 'hadamard', or anything that is listed by 'help wfilters'):
transform_2D_HT_name     = 'bior1.5';%'dct';% %% transform used for the HT filt. of size N1 x N1
transform_2D_Wiener_name = 'dct';     %% transform used for the Wiener filt. of size N1_wiener x N1_wiener
transform_3rd_dim_name   = 'haar';%'dct';%    %% transform used in the 3-rd dim, the same for HT and Wiener filt.

Nf                  = 32;
useM                = 0;
Kin                 = 4;

%%%% Hard-thresholding (HT) parameters:
N1                  = 8;%   %% N1 x N1 is the block size used for the hard-thresholding (HT) filtering
Nstep               = 3;   %% sliding step to process every next reference block
N2                  = 16;  %% maximum number of similar blocks (maximum size of the 3rd dimension of a 3D array)
Ns                  = 39;  %% length of the side of the search neighborhood for full-search block-matching (BM), must be odd
tau_match           = 3000;%% threshold for the block-distance (d-distance)
lambda_thr2D        = 0;   %% threshold parameter for the coarse initial denoising used in the d-distance measure
lambda_thr3D        = 2.7; %% threshold parameter for the hard-thresholding in 3D transform domain
beta                = 2.0; %% parameter of the 2D Kaiser window used in the reconstruction
gamma               = 3.0; %% blockmatching confidence interval
%%%% Wiener filtering parameters:
N1_wiener           = 8;%4;%
Nstep_wiener        = 3;
N2_wiener           = 32;%8;%
Ns_wiener           = 39;
tau_match_wiener    = 400;
beta_wiener         = 2.0;

if strcmp(profile, 'lc') == 1n

    Nstep               = 6;
    Ns                  = 25;
    Nstep_wiener        = 5;
    N2_wiener           = 16;
    Ns_wiener           = 25;
end

% Profile 'vn' was proposed in
%  Y. Hou, C. Zhao, D. Yang, and Y. Cheng, 'Comment on "Image Denoising by Sparse 3D Transform-Domain
%  Collaborative Filtering"', accepted for publication, IEEE Trans. on Image Processing, July, 2010.
% as a better alternative to that initially proposed in [1] (which is currently in profile 'vn_old')
if strcmp(profile, 'vn') == 1

    N2                  = 32;
    Nstep               = 4;

    N1_wiener           = 11;
    Nstep_wiener        = 6;

    lambda_thr3D        = 2.8;
    tau_match_wiener    = 3500;

    Ns_wiener           = 39;

end

% The 'vn_old' profile corresponds to the original parameters for strong noise proposed in [1].
if strcmp(profile, 'vn_old') == 1

    transform_2D_HT_name = 'dct';

    N1                  = 12;
    Nstep               = 4;

    N1_wiener           = 11;
    Nstep_wiener        = 6;

    lambda_thr3D        = 2.8;
    lambda_thr2D        = 2.0;
    tau_match_wiener    = 3500;
    tau_match           = 5000;

    Ns_wiener           = 39;

end

decLevel = 0;        %% dec. levels of the dyadic wavelet 2D transform for blocks (0 means full decomposition, higher values decrease the dec. number)


if strcmp(profile, 'high') == 1 %% this profile is not documented in [1]

    decLevel     = 1;
    Nstep        = 2;
    Nstep_wiener = 2;
    lambda_thr3D = 2.5;
    beta         = 2.5;
    beta_wiener  = 1.5;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Note: touch below this point only if you know what you are doing!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if Nf > 0
   tau_match = 3000;  % Adjust for variance subtraction
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Create transform matrices, etc.
%%%%
[Tfor, Tinv]   = getTransfMatrix(N1, transform_2D_HT_name, decLevel);     %% get (normalized) forward and inverse transform matrices
[TforW, TinvW] = getTransfMatrix(N1_wiener, transform_2D_Wiener_name, 0); %% get (normalized) forward and inverse transform matrices

if ~useM&&((strcmp(transform_3rd_dim_name, 'haar') == 1) || (strcmp(transform_3rd_dim_name(end-2:end), '1.1') == 1))
    %%% If Haar is used in the 3-rd dimension, then a fast internal transform is used, thus no need to generate transform
    %%% matrices.
    hadper_trans_single_den         = {};
    inverse_hadper_trans_single_den = {};
else
    %%% Create transform matrices. The transforms are later applied by
    %%% matrix-vector multiplication for the 1D case.
    for hpow = 0:ceil(log2(max(N2,N2_wiener)))
        h = 2^hpow;
        [Tfor3rd, Tinv3rd]   = getTransfMatrix(h, transform_3rd_dim_name, 0);
        hadper_trans_single_den{h}         = single(Tfor3rd);
        inverse_hadper_trans_single_den{h} = single(Tinv3rd');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 2D Kaiser windows used in the aggregation of block-wise estimates
%%%%
if beta_wiener==2 && beta==2 && N1_wiener==8 && N1==8 % hardcode the window function so that the signal processing toolbox is not needed by default
    Wwin2D = [ 0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.4325    0.6717    0.8644    0.9718    0.9718    0.8644    0.6717    0.4325;
        0.3846    0.5974    0.7688    0.8644    0.8644    0.7688    0.5974    0.3846;
        0.2989    0.4642    0.5974    0.6717    0.6717    0.5974    0.4642    0.2989;
        0.1924    0.2989    0.3846    0.4325    0.4325    0.3846    0.2989    0.1924];
    Wwin2D_wiener = Wwin2D;
else
    Wwin2D           = kaiser(N1, beta) * kaiser(N1, beta)'; % Kaiser window used in the aggregation of the HT part
    Wwin2D_wiener    = kaiser(N1_wiener, beta_wiener) * kaiser(N1_wiener, beta_wiener)'; % Kaiser window used in the aggregation of the Wiener filt. part
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% convert z to double precision if needed
z = double(z);

if numel(PSD)==1
    PSD = ones(size(z))*PSD^2*numel(z);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 1. Produce the basic estimate by HT filtering
%%%%

y_hat = bm3d_thr_colored_noise(z, hadper_trans_single_den, inverse_hadper_trans_single_den, Nstep, N1, N2, ...
        lambda_thr3D, tau_match*N1*N1/(255*255), (Ns-1)/2, single(Tfor), single(Tinv)',...
        zeros(N1, N1), Wwin2D, single(PSD), Nf, gamma, Kin);

disp('Hard-thresholding phase completed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Step 2. Produce the final estimate by Wiener filtering (using the
%%%%  hard-thresholding initial estimate)
%%%

y_est = bm3d_wiener_colored_noise(z, y_hat, hadper_trans_single_den, inverse_hadper_trans_single_den, Nstep_wiener, N1_wiener, N2_wiener, ...
        tau_match_wiener*N1_wiener*N1_wiener/(255*255), (Ns_wiener-1)/2, zeros(N1_wiener, N1_wiener), single(TforW), single(TinvW)',...
        Wwin2D_wiener, single(PSD), Nf, Kin);

disp('Wiener phase completed')

y_est = double(y_est);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some auxiliary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
% Create forward and inverse transform matrices, which allow for perfect
% reconstruction. The forward transform matrix is normalized so that the
% l2-norm of each basis element is 1.
%
% [Tforward, Tinverse] = getTransfMatrix (N, transform_type, dec_levels)
%
%  INPUTS:
%
%   N               --> Size of the transform (for wavelets, must be 2^K)
%
%   transform_type  --> 'dct', 'dst', 'hadamard', or anything that is
%                       listed by 'help wfilters' (bi-orthogonal wavelets)
%                       'DCrand' -- an orthonormal transform with a DC and all
%                       the other basis elements of random nature
%
%   dec_levels      --> If a wavelet transform is generated, this is the
%                       desired decomposition level. Must be in the
%                       range [0, log2(N)-1], where "0" implies
%                       full decomposition.
%
%  OUTPUTS:
%
%   Tforward        --> (N x N) Forward transform matrix
%
%   Tinverse        --> (N x N) Inverse transform matrix
%

if ~exist('dec_levels','var')
    dec_levels = 0;
end

if N == 1
    Tforward = 1;
elseif strcmp(transform_type, 'hadamard') == 1
    Tforward    = hadamard(N);
elseif (N == 8) && strcmp(transform_type, 'bior1.5')==1 % hardcoded transform so that the wavelet toolbox is not needed to generate it
    Tforward =[ 0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110   0.343550200747110
               -0.225454819240296  -0.461645582253923  -0.461645582253923  -0.225454819240296   0.225454819240296   0.461645582253923   0.461645582253923   0.225454819240296
                0.569359398342840   0.402347308162280  -0.402347308162280  -0.569359398342840  -0.083506045090280   0.083506045090280  -0.083506045090280   0.083506045090280
               -0.083506045090280   0.083506045090280  -0.083506045090280   0.083506045090280   0.569359398342840   0.402347308162280  -0.402347308162280  -0.569359398342840
                0.707106781186550  -0.707106781186550                   0                   0                   0                   0                   0                   0
                                0                   0   0.707106781186550  -0.707106781186550                   0                   0                   0                   0
                                0                   0                   0                   0   0.707106781186550  -0.707106781186550                   0                   0
                                0                   0                   0                   0                   0                   0   0.707106781186550  -0.707106781186550];
elseif (N == 8) && strcmp(transform_type, 'dct')==1 % hardcoded transform so that the signal processing toolbox is not needed to generate it
    Tforward = [ 0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274   0.353553390593274;
        0.490392640201615   0.415734806151273   0.277785116509801   0.097545161008064  -0.097545161008064  -0.277785116509801  -0.415734806151273  -0.490392640201615;
        0.461939766255643   0.191341716182545  -0.191341716182545  -0.461939766255643  -0.461939766255643  -0.191341716182545   0.191341716182545   0.461939766255643;
        0.415734806151273  -0.097545161008064  -0.490392640201615  -0.277785116509801   0.277785116509801   0.490392640201615   0.097545161008064  -0.415734806151273;
        0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274   0.353553390593274  -0.353553390593274  -0.353553390593274   0.353553390593274;
        0.277785116509801  -0.490392640201615   0.097545161008064   0.415734806151273  -0.415734806151273  -0.097545161008064   0.490392640201615  -0.277785116509801;
        0.191341716182545  -0.461939766255643   0.461939766255643  -0.191341716182545  -0.191341716182545   0.461939766255643  -0.461939766255643   0.191341716182545;
        0.097545161008064  -0.277785116509801   0.415734806151273  -0.490392640201615   0.490392640201615  -0.415734806151273   0.277785116509801  -0.097545161008064];
elseif (N == 8) && strcmp(transform_type, 'dst')==1 % hardcoded transform so that the PDE toolbox is not needed to generate it
    Tforward = [ 0.161229841765317   0.303012985114696   0.408248290463863   0.464242826880013   0.464242826880013   0.408248290463863   0.303012985114696   0.161229841765317;
        0.303012985114696   0.464242826880013   0.408248290463863   0.161229841765317  -0.161229841765317  -0.408248290463863  -0.464242826880013  -0.303012985114696;
        0.408248290463863   0.408248290463863                   0  -0.408248290463863  -0.408248290463863                   0   0.408248290463863   0.408248290463863;
        0.464242826880013   0.161229841765317  -0.408248290463863  -0.303012985114696   0.303012985114696   0.408248290463863  -0.161229841765317  -0.464242826880013;
        0.464242826880013  -0.161229841765317  -0.408248290463863   0.303012985114696   0.303012985114696  -0.408248290463863  -0.161229841765317   0.464242826880013;
        0.408248290463863  -0.408248290463863                   0   0.408248290463863  -0.408248290463863                   0   0.408248290463863  -0.408248290463863;
        0.303012985114696  -0.464242826880013   0.408248290463863  -0.161229841765317  -0.161229841765317   0.408248290463863  -0.464242826880013   0.303012985114696;
        0.161229841765317  -0.303012985114696   0.408248290463863  -0.464242826880013   0.464242826880013  -0.408248290463863   0.303012985114696  -0.161229841765317];
elseif strcmp(transform_type, 'dct') == 1
    Tforward    = dct(eye(N));
elseif strcmp(transform_type, 'dst') == 1
    Tforward    = dst(eye(N));
elseif strcmp(transform_type, 'DCrand') == 1
    x = randn(N); x(1:end,1) = 1; [Q,~] = qr(x);
    if (Q(1) < 0)
        Q = -Q;
    end;
    Tforward = Q';
else %% a wavelet decomposition supported by 'wavedec'
    %%% Set periodic boundary conditions, to preserve bi-orthogonality
    dwtmode('per','nodisp');

    Tforward = zeros(N,N);
    for i = 1:N
        Tforward(:,i)=wavedec(circshift([1 zeros(1,N-1)],[dec_levels i-1]), log2(N), transform_type);  %% construct transform matrix
    end
end

%%% Normalize the basis elements
if ~((N == 8) && strcmp(transform_type, 'bior1.5')==1)
    Tforward = (Tforward' * diag(sqrt(1./sum(Tforward.^2,2))))';
end

%%% Compute the inverse transform matrix
Tinverse = inv(Tforward);

return;

