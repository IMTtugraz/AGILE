function [ v ] = readbin_vector( filename )
%READBIN_CRSMATRIX reads a sparse matrix from a binary crs matrix file
%   If filename is not set the file can be choosen with an open dialog

    v = [];
    path = '';

    if nargin == 0
        [file, path] = uigetfile({'*.bin','Binary CRS matrix file (*.bin)';'*.*',  'All Files (*.*)'},'Load crs matrix ...');
        filename = fullfile(path, file);
    end
    
    if isequal(filename,0) || isequal(path,0)
        % user pressed cancel
        return;
    end
    
    % try to open a file for reading (binary)
    [fid, message] = fopen(filename, 'r');
    if fid == -1
        error(message);
    end
    
    % read complex header info
    iscomplex = fread(fid, 1, 'uchar');

    % read matrix size
    rows = fread(fid, 1, 'uint32');   % read number of rows
    cols = fread(fid, 1, 'uint32');   % read number of columns
    nnz = fread(fid, 1, 'uint32');   % read total number of nonzeros

    % read number of nonzeros per row (rows elements)
    row_nnz = fread(fid, rows, 'uint32');
    
    % read column indices (zero based)
    column_ind = fread(fid, nnz, 'uint32') + 1;
    
    % read column indices (zero based)
    data = fread(fid, nnz, 'double');
    if bitand(iscomplex, 1) > 0
        data = complex(data, fread(fid, nnz, 'double'));
    end

    row_ind = zeros(nnz,1);
    row_nnz = cumsum([1;row_nnz]);
    for ii = 1:rows
        row_ind(row_nnz(ii):row_nnz(ii+1)-1) = ii;
    end
    
    v = sparse(row_ind, column_ind, data, rows, cols, nnz);
    
    fclose(fid);
    
end
