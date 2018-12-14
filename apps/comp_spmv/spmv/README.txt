            IBM Sparse Matrix-Vector Multiplication Toolkit for Graphics Processing Units

Requirements for installing SpMV Toolkit
=========================================

1. A x86 or x86-64 CPU which is connected to a "CUDA-enabled" GPU device
(e.g., GeForce 8800 GTX, 8800 GTS, GTX 280, GTX 260).
2. Linux OS only
3. Machine should support 32-bit ELF file format.
4. CUDA driver and CUDA toolkit should have been installed. Please visit
http://www.nvidia.com/object/cuda_get.html for more details on how to download and install the CUDA driver and toolkit for your machine.
5. The toolkit has been tested on CUDA 2.0 and 2.1 on 8800 GTX and GTX 280 GPUs.

Files contained in the SpMV Toolkit
===================================

README.txt (this file)
FAQ
licences sub-directory
paper sub-directory containing IBM Research TR 24704
Makefile
cache.h config.h SpMV.h SpMV_compute.cpp SpMV_inspect.h   
SpMV_alloc.cu SpMV_gen.cpp SpMV_inspect.cpp  comp_SpMV.cu  SpMV_compute.cpp 

Installing SpMV Toolkit 
=======================

1. Download the package 'SpMV-src-package.tar.gz'. 

2. Untar the package and invoke make
	$ tar zxvf SpMV-package.tar.gz
        $ make

3. Let the CUDA install path on the machine be <CUDA_INSTALL_PATH>. Include
<CUDA_INSTALL_PATH>/lib in LD_LIBRARY_PATH environment varibale, so that the
executable can locate the 'cuda runtime library' while execution.  

4. Run the executable with the following command line arguments
	a.'input sparse matrix' path - mandatory
	b.'input vector' path - optional
	c.'output vector' path - optional [if given, 'input vector' path
should also be given]
	
	$ ./SpMV <input matrix path> [ <input vector path> [ <output vector
path> ] ]


For further information, please contact Rajesh Bordawekar at
bordaw@us.ibm.com.
