% Skript: Test 3d Array File I/O with IMT CUDA
% Last Change: 3.2.2010, 15:00, by Florian

addpath('/home/knoll/svn/imtcuda/include/imtcuda/io');
addpath('/home/knoll/Documents/MATLAB/Optimized_Random_IRGN/utils');

clear all; close all; clc;

showCoil = 4;

%% load data
load t2brainRawdata.mat

img = ifft2c(rawdata);
img = single(img);

[num_rows, num_columns, num_coils] = size(img);
num_bytes_per_entry = 4;
img_is_complex = ~isreal(img);

%% write image data to a file
writeMatlab2Bin_3D(single(img), 'img.bin');

%% perform fft on GPU with cuda
% system('./3D_Matrix_I_O img.bin img2.bin');

%% read in the result
cuda_recon = readBin2Matlab_3D('img2.bin', num_rows, num_columns, num_coils, num_bytes_per_entry);

%% Compare the results: do some plots
figure, 
subplot(1,2,1), imshow(abs(img(:,:,showCoil)),[]);
subplot(1,2,2), imshow(abs(cuda_recon(:,:,showCoil)),[]);
% subplot(1,3,3), imshow(abs(img(:,:,showCoil)-cuda_recon(:,:,showCoil)),[]);
