function writeMatlab2Bin_3D(matrix, filename)
% WRITEMATLAB2BIN_3D(MATRIX,FILENAME)
% Writes a Matlab matrix to a bin file that can be accessed with AGILE
%
% INPUT:
% matrix: any 3D matrix, Single precision, real or complex
% filename: string
%
% OUTPUT:
% None, but a bin file is written in the current directory
%
% Example: writeMatlab2Bin_3D(testmatrix, 'testfile.bin');
%
% See also: readBin2Matlab_3D.m
%
% Last Change: 4.2.2010, 13:00
% By: Florian


%% Prepare some info variables
matrix = single(matrix);
[numRows, numColumns, numCoils] = size(matrix);
numBytesPerEntry = 4;
matrixIsComplex = ~isreal(matrix);

%% Start output prodcedure
fid = fopen(filename, 'wb');
if fid == -1
    error(['Could not open file ', filename '']);
end

% write real matrix to a file
if ~matrixIsComplex
    fwrite(fid, numRows, 'uint32');
    fwrite(fid, numColumns, 'uint32');
    fwrite(fid, numCoils, 'uint32');
    fwrite(fid, numBytesPerEntry, 'uint32');
    fwrite(fid, matrixIsComplex, 'uint32');
    for ii=1:numCoils
        fwrite(fid, matrix(:,:,ii).', 'float32');
    end
    fclose(fid);
else
    % write complex matrix to a file
    fwrite(fid, numRows, 'uint32');
    fwrite(fid, numColumns, 'uint32');
    fwrite(fid, numCoils, 'uint32');
    fwrite(fid, numBytesPerEntry, 'uint32');
    fwrite(fid, matrixIsComplex, 'uint32');
    
    % Seperate real and imaginary part
    tempMatrix = zeros(numRows,numColumns*2,numCoils);
    tempMatrix(:, 1:2:end-1,:) = real(matrix);
    tempMatrix(:, 2:2:end,:) = imag(matrix);
    for ii=1:numCoils
        fwrite(fid, tempMatrix(:,:,ii).', 'float32');
    end
    fclose(fid);    
end
        
        
