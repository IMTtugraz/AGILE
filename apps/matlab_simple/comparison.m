num_rows = 2000;
num_columns = 2000;
num_bytes_per_entry = 4;
A_is_complex = 0;
x_is_complex = 0;

% create a dense matrix
A = (rand(num_rows, num_columns) - 0.5) * 1000;
if A_is_complex
  A = A + i * (rand(num_rows, num_columns) - 0.5) * 1000;
end
A = single(A);
% create a vector
x = (rand(num_columns, 1) - 0.5) * 1000;
if x_is_complex
  x = x + i * (rand(num_columns, 1) - 0.5) * 1000;
end
x = single(x);

tic;
y = A * x;
matlab_time = toc;

% write the matrix to a file
fid = fopen('A.bin', 'wb');
if fid == -1
  error('Could not open file ''A.bin''');
end
fwrite(fid, num_rows, 'uint32');
fwrite(fid, num_columns, 'uint32');
fwrite(fid, num_bytes_per_entry, 'uint32');
fwrite(fid, A_is_complex, 'uint32');
fwrite(fid, A.', 'float32');
fclose(fid);

% write x to a file
fid = fopen('x.bin', 'wb');
if fid == -1
  error('Could not open file ''x.bin''');
end
fwrite(fid, num_columns, 'uint32');
fwrite(fid, num_bytes_per_entry, 'uint32');
fwrite(fid, x_is_complex, 'uint32');
fwrite(fid, x, 'float32');
fclose(fid);

% perform the multiplication
tic
system('./dense_multiplication.exe A.bin x.bin y.bin');
cuda_time = toc;

% read in the result
fid = fopen('y.bin', 'rb');
if fid == -1
  error('Could not open file ''y.bin''');
end
y_info.size = fread(fid, 1, 'uint32');
if y_info.size ~= num_rows
  error('Dimension mismatch');
end
y_info.num_byte_per_entry = fread(fid, 1, 'uint32');
if y_info.num_byte_per_entry ~= num_bytes_per_entry
  error('Type mismatch');
end
y_info.is_complex = fread(fid, 1, 'uint32');
if y_info.is_complex
  temp = fread(fid, y_info.size * 2, 'float32');
  cuda_y = temp(1 : 2 : end) + i * temp(1 : 2 : end);
else
  cuda_y = fread(fid, y_info.size, 'float32');
end
fclose(fid);

% compare the result
abs_difference = abs(y - cuda_y);
indices = y ~= 0;
absolute_error = max(abs_difference);
relative_error = max(abs_difference(indices) ./ y(indices));

disp(['Matlab time: ', num2str(matlab_time)]);
disp(['CUDA time:   ', num2str(cuda_time)]);
disp(['Max absolute error: ', num2str(absolute_error)]);
disp(['Max relative error: ', num2str(relative_error * 100), '%']);
disp(' ');
disp(' ');
