INCLUDES=-I.
LIBS=-lm -lfftw3f -lfftw3f_threads -lpthread
LIBS_NO_NFFT=-lm -lfftw3f -lfftw3f_threads -lpthread

CC=gcc
CXX=g++
CFLAGS=-g -O3 -fopenmp -ffast-math $(INCLUDES) `pkg-config gtk+-2.0 --cflags`
LINKFLAGS=-fopenmp $(LIBPATH)

default: nfft2d_tgv_recon_CPU

nfft2d_tgv_recon_CPU: nfft2d_tgv_recon_CPU.o nfft2d.o ParseParam.o
	$(CXX) -o nfft2d_tgv_recon_CPU nfft2d_tgv_recon_CPU.o nfft2d.o ParseParam.o $(LINKFLAGS) $(LIBS) 

nfft2d_tgv_recon_CPU.o: nfft2d_tgv_recon_CPU.cpp
	$(CXX) $(CFLAGS) -c nfft2d_tgv_recon_CPU.cpp

nfft2d.o: nfft2d.c
	$(CXX) -c $(CFLAGS) nfft2d.c

ParseParam.o: ParseParam.c
	$(CXX) $(CFLAGS) -c ParseParam.c

clean:
	rm nfft2d.o ParseParam.o nfft2d_tgv_recon_CPU.o nfft2d_tgv_recon_CPU

test: nfft2d_tgv_recon_CPU
	./nfft2d_tgv_recon_CPU --param recon_24_param.txt
