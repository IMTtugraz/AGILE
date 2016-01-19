function [ ] = writebin_crsmatrix( A, filename )
%WRITEBIN_CRSMATRIX Saves a matrix in binary crs format
%   If filename is not set a save dialog will be opened

    path = '';

    if nargin < 1
        error('Not enough input arguments.');
    elseif nargin == 1
        [file, path] = uiputfile({'*.bin','Binary CRS matrix file (*.bin)';'*.*',  'All Files (*.*)'},'Save crs matrix ...');
        filename = fullfile(path, file);
    end

    if isequal(filename,0) || isequal(path,0)
        % user pressed cancel
        return;
    end
    
    % try to open a file for writing (binary)
    [fid, message] = fopen(filename, 'w');
    if fid == -1
        error(message);
    end

    % ensure double type
    A = double(A);
    
    % check for complex numbers
    iscomplex = 0;
    if (~isreal(A))
        iscomplex = 1;
    end
    
    % write file header information for complex
    fwrite(fid, iscomplex, 'uchar');
    
    % compute data
    [rows, cols] = size(A);
    row_nnz = full(sum(A~=0,2));  % number of non zeros per row
    [column_index, row_index, data] = find(A.');

    % start writing data
    fwrite(fid, rows, 'uint32');          % write number of rows
    fwrite(fid, cols, 'uint32');          % write number of columns
    fwrite(fid, size(data,1), 'uint32');  % write total number of non zeros
       
    % write number of non zeros per row (rows elements)
    fwrite(fid, row_nnz, 'uint32');

    %write column index (zero based)
    fwrite(fid, (column_index-1), 'uint32');
    
    % write data coresponding to the column index
    fwrite(fid, real(data), 'double');
    if iscomplex
        % negative sign because of A' above!
        fwrite(fid, imag(data), 'double');
    end
    
    fclose(fid);
    
    disp(sprintf('Matrix "%s" saved to: "%s"', inputname(1), filename));
end


