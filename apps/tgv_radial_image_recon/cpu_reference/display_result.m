% Generate plot of CPU reference
clear all; close all; clc;

% Display Parameters
windowing = 0.8;
brightening = 0.2;

% 24 projections load TGV2 result
fid = fopen('./recon_24_result.dat');
recon = fread(fid,inf, 'float32');
fclose(fid);
n = size(recon,1);
TGV2_recon_24proj = reshape(recon,[sqrt(n) sqrt(n) n/sqrt(n)^2]);
TGV2_recon_24proj = abs(TGV2_recon_24proj)';
TGV2_recon_24proj = TGV2_recon_24proj/max(TGV2_recon_24proj(:));

% Plot
figure, imshow(brighten(TGV2_recon_24proj,brightening),[0, windowing]);
title('TGV CPU Reference');