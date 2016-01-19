% Construct radial rawdata for CUDA, starting from a 96 radial projections
% rawdata set

clear all; close all; clc;

%% Load Data
load rawdata_96proj_4ch.mat;
dataRadial = rawdata; clear rawdata;

%% Generate CUDA style data files
for skip=[1,2,3,4]
    
    rawdata = dataRadial;
    rawdata = rawdata(:,skip:skip:end,:);

    datalen = size(rawdata,1);
    spokes = size(rawdata,2);
    coils = size(rawdata,3);

    % compute coordinates
    angles = (1:spokes)/spokes*pi;
    line = linspace(-0.5,0.5,datalen);
    [arg1,arg2] = meshgrid(angles,line);
    coordx = arg2.*cos(arg1);
    coordy = arg2.*sin(arg1);
    coord = [coordx(:)'; coordy(:)'];

    % compute weights
    weights = sqrt(abs(arg2(:))/spokes);

    % compute data
    datareal = real(rawdata(:));
    dataimag = imag(rawdata(:));
    data = [datareal'; dataimag'];

    % write everything
    fid = fopen(sprintf('recon_%02d.crd',spokes),'w');
    fwrite(fid, coord', 'float32');
    fclose(fid);

    fid = fopen(sprintf('recon_%02d.wgt',spokes),'w');
    fwrite(fid, weights, 'float32');
    fclose(fid);

    fid = fopen(sprintf('recon_%02d.dat',spokes),'w');
    fwrite(fid, data, 'float32');
    fclose(fid);
end




