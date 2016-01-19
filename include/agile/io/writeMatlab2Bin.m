function writeMatlab2Bin(matrix, filename)
% WRITEMATLAB2BIN(MATRIX,FILENAME)
% Writes a matlab matrix to a bin file that can be accessed with AGILE
%
% INPUT:
% matrix: any 2D matrix, Single precision, real or complex
% filename: string
%
% OUTPUT:
% None, but a bin file is written in the current directory
%
% Example: writeMatlab2Bin(testmatrix, 'testfile.bin');
%
% See also: readBin2Matlab.m
%
% Last Change: 26.3.2009, 14:00
% By: Florian


%% Prepare some info variables
matrix = single(matrix);
[numRows, numColumns] = size(matrix);
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
    fwrite(fid, numBytesPerEntry, 'uint32');
    fwrite(fid, matrixIsComplex, 'uint32');
    fwrite(fid, matrix.', 'float32');
    fclose(fid);
else
    % write complex matrix to a file
    fwrite(fid, numRows, 'uint32');
    fwrite(fid, numColumns, 'uint32');
    fwrite(fid, numBytesPerEntry, 'uint32');
    fwrite(fid, matrixIsComplex, 'uint32');
    
    % Seperate real and imaginary part
    tempMatrix = zeros(numRows,numColumns*2);
    tempMatrix(:, 1:2:end-1) = real(matrix);
    tempMatrix(:, 2:2:end) = imag(matrix);
    fwrite(fid, tempMatrix.', 'float32');
    fclose(fid);    
end
        
        
