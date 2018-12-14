% Skript: Test performance of CUFFT with brain.mat file
% - Writes a binary file from an MRI brain mat file
% - Runs fft_test.exe which does fft and ifft
% - Reads in result of CUFFT computation and shows matlab and CUFFT results
% 
% NOTE: Calling fft_test.exe might not work due to problems with libstdc++:
% In this case, you have to run fft_test.exe from the command line:
% ./fft_test.exe img.bin img_recon.bin
%
% Last Change: 26.3.2009, 16:00, by Florian

addpath('/home/knoll/svn/imtcuda/include/imtcuda/io'); % Change this

clear all; close all; clc;

%% load data
load img_brain_t2.mat
minim = min(img_brain_t2(:));
maxim = max(img_brain_t2(:));
img = (img_brain_t2-minim)./(maxim-minim);
clear img_brain_t2;
img = img + i*0.1*randn(size(img)); % Just add some random complex stuff
% img = imresize(img,0.5);
img = single(img);

%% Do Padding in order to make matrix size odd
padding = zeros(size(img,1),1);
img = [img, padding];
padding = zeros(1, size(img,2));
img = [padding; img];

[num_rows, num_columns] = size(img);
num_bytes_per_entry = 4;
img_is_complex = ~isreal(img);

%% perform reference 2d fft in matlab
matlab_kspace = fft2(img);
matlab_img_recon = ifft2(matlab_kspace);

%% write image data to a file
writeMatlab2Bin(img, 'img.bin');

%% perform fft on GPU with cuda
% system('./fft_test.exe img.bin img_recon.bin');

%% read in the result
cuda_kspace = readBin2Matlab('kspace.bin', num_rows, num_columns, num_bytes_per_entry);

%% Compare the results: do some plots
figure, imshow(abs(matlab_img_recon),[]);
figure, imshow(log(1+abs(cuda_kspace)),[]);
