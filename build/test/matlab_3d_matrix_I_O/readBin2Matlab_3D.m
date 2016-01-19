function matrix = readBin2Matlab_3D(filename, numRows, numColumns, numCoils, numBytesPerEntry)
% MATRIX = READBIN2MATLAB_3D(FILENAME, NUMROWS, NUMCOLUMS, NUMBYTESPERENTRY)
% Reads a bin file from AGILE in Matlab
%
% INPUT:
% filename: string, name of the bin file
% numRows 
% numColums
% numColums
% numBytesPerEntry
%
% OUTPUT:
% matrix: any 3D matrix, Single precision, real or complex
%
% Example: testMatrix = readBin2Matlab_3D('testfile1.bin', numRows,
%                       numColums, numBytesPerEntry);
%
% See also: writeMatlab2bin_3D.m
%
% Last Change: 4.2.2010, 13:00
% By: Florian

%% open file
fid = fopen(filename, 'rb');
if fid == -1
  error(['Could not open file ', filename '']);
end

%% Check file information
matrix_info.rows = fread(fid, 1, 'uint32');
if matrix_info.rows ~= numRows
  error('Dimension mismatch');
end
matrix_info.columns = fread(fid, 1, 'uint32');
if matrix_info.columns ~= numColumns
  error('Dimension mismatch');
end
matrix_info.coils = fread(fid, 1, 'uint32');
if matrix_info.coils ~= numCoils
  error('Dimension mismatch');
end
matrix_info.num_byte_per_entry = fread(fid, 1, 'uint32');
if matrix_info.num_byte_per_entry ~= numBytesPerEntry
  error('Type mismatch');
end

%% Determine if file contains complex or real content
matrix_info.is_complex = fread(fid, 1, 'uint32');
if matrix_info.is_complex
  temp = single(fread(fid, matrix_info.rows * matrix_info.columns * 2 * matrix_info.coils, 'float32'));
  matrix = temp(1 : 2 : end) + i * temp(2 : 2 : end);
  clear temp;
else
  matrix = single(fread(fid, matrix_info.rows * matrix_info.columns * matrix_info.coils, 'float32'));
end

%% Arrange in the right shape: Note that matlab permutes rows and columns
%% in comparison to C. This is why the reshape command has to be called that way 
fclose(fid);
matrix = reshape(matrix,numColumns,numRows,numCoils);
for ii = 1:numCoils
    temp(:,:,ii) = matrix(:,:,ii).';
end
matrix = temp; clear temp;
