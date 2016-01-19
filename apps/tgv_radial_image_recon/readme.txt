This application performs second order total generalized variation (TGV2) image reconstruction from multiple coil data, as described in [1].

You need a k-space data file (e.g. recon_32.dat), a file describing the k-space trajectory (e.g. recon_32.crd) and if you want to use any, a file with density compensation weights (e.g. recon_32.wgt). The best way is probably to take a look at the file make_CUDA_data_radial_4ch.m to see how these are organised.

You also need a file defining the reconstruction parameters (e.g. recon_32_param.txt). This should be pretty straight forward to adapt to your data.

After compiling the library, the executeable is found in the ../build/apps/tgv_radial_image_recon (wherever you created your build directory). The reconstruction is then started with:

./nfft2d_tgv_recon_nogfx -param=your_parameter_file.txt (e.g. ./nfft2d_tgv_recon_nogfx -param=recon_32_param.txt)

The corresponding data and param files have to be in the same directory as the executable. If you want to test the code with the provided test data, just copy the provided *.dat, *.crd, *.wgt and *param.txt files to the build directory.

A version of the reconstruction code with online display of the current iteration is also available. However, at the moment, this version is not included in the standard version of the library as it requires gtkimageviewer, and this poses some problems with compilation when integrated in the rest of the library.

The folder "cpu_reference" contains a CPU reference implentation. However, as this is not really considered to be a part of the library, it will not be compiled automatically with cmake. A seperate makefile is included in this folder. To test it, please copy the included data in the "cpu_reference" folder and run with "./nfft2d_tgv_recon_CPU --param recon_32_param.txt" in the case of the 32 projections data.

[1] Knoll, F.; Bredies, K.; Pock, T.; Stollberger, R.: 
Second Order Total Generalized Variation (TGV) for MRI: Magnetic Resonance in Medicine, 65 (2011), P. 480 - 491
