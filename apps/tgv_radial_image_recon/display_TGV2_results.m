% Display CUDA TGV2 results

clear all; close all; clc;

%addpath(genpath('/home/knoll/Documents/MATLAB/NUFFT'));
%addpath('/home/knoll/Documents/MATLAB/Optimized_Random_IRGN/utils');

%% Load TGV2 result
fid = fopen('./recon_48_result.dat');
recon = fread(fid,inf, 'float32');
fclose(fid);
n = size(recon,1);
TGV2_recon = reshape(recon,[sqrt(n) sqrt(n) n/sqrt(n)^2]);
TGV2_recon = TGV2_recon';
nPE = size(TGV2_recon,1);
nFE = size(TGV2_recon,2);
figure, imshow(TGV2_recon,[]), title('TGV2 Reconstruction');

%% Comparison: NUFFT reconstruction

% Load Data
load rawdata_96proj_4ch.mat;
nCh = size(rawdata,3);

skip=3;
rawdata = rawdata(:,skip:skip:end,:);

% Generate NUFFT Operator
spokes = size(rawdata,2);
datalen = size(rawdata,1);

angles = (1:spokes)/spokes*pi;
line = linspace(-0.5,0.5,datalen);
[arg1,arg2] = meshgrid(angles,line);
coordx = arg2.*cos(arg1);
coordy = arg2.*sin(arg1);
k = coordx+i*coordy;

% generate Fourier sampling operator
FT = NUFFT(k, 1, 1, 0, [nPE,nFE], 2);

% Cheap and Simple Version of Density Compensation
w = (abs(line)*datalen/spokes/2)';
w = repmat(w, [1, spokes, nCh]);

% IFT and SOS combination
NUFFT_recon = FT'*(rawdata.*w);
NUFFT_recon  = sqrt(sum(abs(NUFFT_recon ).^2,3));

%% Display results
figure, imshow(NUFFT_recon,[]), title('Min. Norm Reconstruction');
figure, imshow(TGV2_recon,[]), title('TGV2 Reconstruction');
