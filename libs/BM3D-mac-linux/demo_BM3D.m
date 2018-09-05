% close all
% clearvars

y=im2double(imread('cameraman256.png'));
% y=im2double(imread('montage.png'));
% y=im2double(imread('Lena512.png'));

k=eye(8);k=conv2(k,ones(2))/140;
%k=ones(1, 28);k=conv2(k,ones(2))/140;
%k=ones(3,3);k=conv2(k,ones(2))/140;
% k=35/255;

randn('seed',0);

z=y+imfilter(randn(256),k(end:-1:1,end:-1:1),'circular');
PSD=abs(fft2(k,size(z,1),size(z,2))).^2*numel(z);

%%
disp('*------------------------*')
tic
useM=0;
[y_est8, y_est_HT8] = BM3D(z, PSD, 'np');
toc
imshow(y_est8)
disp('*------------------------*')