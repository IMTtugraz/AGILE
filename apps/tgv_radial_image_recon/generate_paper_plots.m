%% Generate TGV2 radial recon plots for IMT CUDA Paper

clear all; close all; clc;

% Display Parameters
windowing = 0.8;
brightening = 0.2;

%% 96 projections load TGV2 result
fid = fopen('./recon_96_result.dat');
recon = fread(fid,inf, 'float32');
fclose(fid);
n = size(recon,1);
TGV2_recon_96proj = reshape(recon,[sqrt(n) sqrt(n) n/sqrt(n)^2]);
TGV2_recon_96proj = abs(TGV2_recon_96proj)';
TGV2_recon_96proj = TGV2_recon_96proj/max(TGV2_recon_96proj(:));

%% 48 projections load TGV2 result
fid = fopen('./recon_48_result.dat');
recon = fread(fid,inf, 'float32');
fclose(fid);
n = size(recon,1);
TGV2_recon_48proj = reshape(recon,[sqrt(n) sqrt(n) n/sqrt(n)^2]);
TGV2_recon_48proj = abs(TGV2_recon_48proj)';
TGV2_recon_48proj = TGV2_recon_48proj/max(TGV2_recon_48proj(:));

%% 32 projections load TGV2 result
fid = fopen('./recon_32_result.dat');
recon = fread(fid,inf, 'float32');
fclose(fid);
n = size(recon,1);
TGV2_recon_32proj = reshape(recon,[sqrt(n) sqrt(n) n/sqrt(n)^2]);
TGV2_recon_32proj = abs(TGV2_recon_32proj)';
TGV2_recon_32proj = TGV2_recon_32proj/max(TGV2_recon_32proj(:));

%% 24 projections load TGV2 result
fid = fopen('./recon_24_result.dat');
recon = fread(fid,inf, 'float32');
fclose(fid);
n = size(recon,1);
TGV2_recon_24proj = reshape(recon,[sqrt(n) sqrt(n) n/sqrt(n)^2]);
TGV2_recon_24proj = abs(TGV2_recon_24proj)';
TGV2_recon_24proj = TGV2_recon_24proj/max(TGV2_recon_24proj(:));


%% Comparison: NUFFT reconstruction
% Load Data
load rawdata_96proj_4ch.mat;
rawdata_96proj = rawdata; clear rawdata;
nCh = size(rawdata_96proj,3);
nPE = size(TGV2_recon_24proj,1);
nFE = size(TGV2_recon_24proj,2);

NUFFT_recon_sos = zeros(nPE,nFE,nCh);
for skip = [1,2,3,4]
    
    rawdata = rawdata_96proj(:,skip:skip:end,:);

    % Generate NUFFT Operator
    spokes = size(rawdata,2);
    datalen = size(rawdata,1);

    angles = (1:spokes)/spokes*pi;
    line = linspace(-0.5,0.5,datalen);
    [arg1,arg2] = meshgrid(angles,line);
    coordx = arg2.*cos(arg1);
    coordy = arg2.*sin(arg1);
    k = coordx+1i*coordy;

    % generate Fourier sampling operator
    FT = NUFFT(k, 1, 1, 0, [nPE,nFE], 2);

    % Cheap and Simple Version of Density Compensation
    w = (abs(line)*datalen/spokes/2)';
    w = repmat(w, [1, spokes, nCh]);

    % IFT and SOS combination
    NUFFT_recon = FT'*(rawdata.*w);
    NUFFT_recon = abs(sqrt(sum(abs(NUFFT_recon ).^2,3)));
    NUFFT_recon_sos(:,:,skip) = NUFFT_recon/max(NUFFT_recon(:));
end

%% Display results
spacing = ones(size(NUFFT_recon_sos,1),1);

results_96 = [NUFFT_recon_sos(:,:,1), spacing, TGV2_recon_96proj];
results_48 = [NUFFT_recon_sos(:,:,2), spacing, TGV2_recon_48proj];
results_32 = [NUFFT_recon_sos(:,:,3), spacing, TGV2_recon_32proj];
results_24 = [NUFFT_recon_sos(:,:,4), spacing, TGV2_recon_24proj];

spacing2 = ones(40, size(results_96,2));
spacing3 = ones(1, size(results_96,2));

results = [spacing2; results_96; spacing3; results_48; spacing3; results_32; spacing3; results_24];
figure, imshow(brighten(results,brightening),[0, windowing]);

%% Add text
h = text(size(TGV2_recon_24proj,2)/2,22,'Min. Norm');
set(h,'FontSize',14,'Horizontalalignment','center');
h = text(size(TGV2_recon_24proj,2)*3/2+2,22,'TGV');
set(h,'FontSize',14,'Horizontalalignment','center');

h = text(size(TGV2_recon_24proj,2)*2-20,20+size(TGV2_recon_24proj,2),'96');
set(h,'FontSize',14,'Horizontalalignment','center','VerticalAlignment','middle','Color', [1 1 1]);
h = text(size(TGV2_recon_24proj,2)*2-20,20+size(TGV2_recon_24proj,2)*2,'48');
set(h,'FontSize',14,'Horizontalalignment','center','VerticalAlignment','middle','Color', [1 1 1]);
h = text(size(TGV2_recon_24proj,2)*2-20,20+size(TGV2_recon_24proj,2)*3,'32');
set(h,'FontSize',14,'Horizontalalignment','center','VerticalAlignment','middle','Color', [1 1 1]);
h = text(size(TGV2_recon_24proj,2)*2-20,20+size(TGV2_recon_24proj,2)*4,'24');
set(h,'FontSize',14,'Horizontalalignment','center','VerticalAlignment','middle','Color', [1 1 1]);
