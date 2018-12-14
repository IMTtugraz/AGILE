# AGILE
Environment for Linear and non-linear Image reconstruction using GPU Acceleration

## Contributors (sorted alphabetically):
* [Kristian Bredies](http://www.uni-graz.at/~bredies) (University of Graz) 
* Gerald Buchgraber (Graz University of Technology, now at Datenkraft IT-Consulting)
* [Manuel Freiberger](mailto:manuel.freiberger@gmx.at) (Graz University of Technology, now at Anton Paar)
* Andreas Huber (Graz University of Technology)
* [Florian Knoll](http://cai2r.net/people/florian-knoll) (CAI2R, New York School of Medicine)


## General Information
AGILE (Environment for Linear and non-linear Image reconstruction using Gpu
Acceleration) is an open source library for GPU accelerated reconstruction
problems in medical imaging. It includes reconstruction code for Fluorescence
Tomography and Magnetic Resonance Imaging.

If you use this library or the provided data in your publications, please cite:
Knoll, F.; Freiberger, M.; Bredies, K.; Stollberger, R.: AGILE: An open source library
for image reconstruction using graphics card hardware acceleration. Proc. Intl. Soc. Mag.
Reson. Med. 19:2554 (2011)

This work is funded and supported by the Austrian Science Fund (FWF) in the context of project 'SFB F3209-19' (Mathematical Optimization and Applications in Biomedical Sciences)
[MOBIS](http://math.uni-graz.at/mobis/)

## Build & Installation
AGILE was developed and tested on Linux (OpenSUSE and Ubuntu) systems. It uses
CMake in order to generate the Makefiles in a flexible manner. You also need
the third-party libraries Boost (http://www.boost.org/) and OpenMPI
(http://www.open-mpi.org/) to compile the library.

It is suggested to make an out-of-source build where the resultant object
files, libraries and executables will be placed separated from the source
files. This helps to keep the source directories clean.

To build the library, create a new directory (e.g. 'build') and execute the
'cmake' command in there providing the root directory of AGILE (which should
be the current directory) as argument. Afterwards 'make' can be executed
in the build directory to generate the library.

Copy/paste instructions for the lazy to be executed right in this directory:
``` 
 mkdir build
 cd build
 cmake ..
``` 

At this step you need to set the path where you installed the CUDA SDK:
``` 
 ccmake .
``` 

Toggle to advanced mode (press 't') and enter the path where you installed the
CUDA SDK in the variable "CUDA_SDK_ROOT_DIR". At this point, you can also
specify to enable double-precision floating point support on the GPU. This
can be done by setting the variable ENABLE_GPU_DOUBLE to ON. Note that this
is a requirement to build the TGV-MRI example which will be skipped if the flag
is not set.

Then compile the library using
``` 
make
``` 

####Note:
AGILE was tested on 64 bit platforms. If you want to compile it on a
32 bit system, you have to change two entries in
/apps/tgv_radial_image/recon/CMakeLists.txt
"link_libraries("-L/${CUDA_SDK_ROOT_DIR}/C/lib -lcutil_x86_64")" to
"link_libraries("-L/${CUDA_SDK_ROOT_DIR}/C/lib -lcutil")"
"target_link_libraries(nfft2d_tgv_recon_nogfx ${CUDA_cufft_LIBRARY} cutil_x86_64)"
to "target_link_libraries(nfft2d_tgv_recon_nogfx ${CUDA_cufft_LIBRARY} cutil)"

If CMake cannot determine the location of the CUDA toolkit, you have to
provide an variable called CUDA_TOOLKIT_ROOT_DIR that is set to the root
directory of the CUDA toolkit. Normally this is the parent directory
of the directory in which nvcc, the CUDA compiler, resides in. The following
commands can be used on a system where the full path to the CUDA compiler
is /opt/cuda/bin/nvcc:
``` 
mkdir build
cd build
cmake -DCUDA_TOOLKIT_ROOT_DIR=/opt/cuda ..
make
``` 


To install the compiled library, simply call
``` 
make install
``` 


On Linux machines the directory defaults to /usr/local/.

You can specify a different directory for installation by compiling and
installing the library with
``` 
cmake -DCMAKE_INSTALL_PREFIX=/prefix/for/the/library ..
make install
``` 


The commands above will install the files into the directory
/prefix/for/the/library/.

## First steps
If you are new to the AGILE library, the tutorials are a good place to start.
Run 'make doc' in the build directory. Afterwards open the file
'doc/html/index.html' in your favourite browser and proceed to the tutorials.

If you are especially interested in iterative TGV based MR image
reconstruction, please see Florian Knolls [AGILE-TGV](http://www.cai2r.net/resources/software/agile-gpu-image-reconstruction-library) implementation.

IRGN-T(G)V for Cartesian image reconstruction:
The IRGN executable is available in the 'bin' folder.
Type './IMT_IRGN --help' to see all available command line options.
Example:
./IMT_IRGN --input=braindata.bin --output=recon.dcm --format=2 --calc=1


## Version history
* version 1.4: 06.06.2018: Added command line functionality for imt_irgn application
* version 1.3: 01.01.2016: Included necessary extensions for dynamic MRI reconstrucion ([AVIONIC](https://github.com/IMTtugraz/AVIONIC)) by [Andreas Schwarzl](https://github.com/andyschwarzl)
  * Refactoring to support latest CUDA architecture
  * Finite forward and backward differences for 3-D data vectors
  * Forward and backward FFT operations for vector based data sets
  * Extended vector utilities (l1/l2-norm, min-max element computation, element-wise compare)
  * Extended thread assignment for large data sets
* version 1.2: 01.01.2014: Included IRGN-T(G)V for Cartesian image reconstruction (apps/imt_irgn) by [Herbert Heigl](https://github.com/hheigl)
* version 1.1: 16.01.2012: Included CPU reference implementation for TGV.
* version 1.0: 28.07.2011: First release of AGILE.
