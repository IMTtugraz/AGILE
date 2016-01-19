
CUDA_BASE_PATH = /usr/local/cuda/
#CUDA_BASE_PATH = /opt/cuda/

LIB_PATHS = -L$(CUDA_BASE_PATH)lib
INCL_PATHS = -I$(CUDA_BASE_PATH)include
OPTFLAGS = -O3
CC = g++
#NVCC = nvcc -ccbin ./
NVCC = nvcc 
OBJS = comp_SpMV.o SpMV_alloc.o SpMV_compute.o SpMV_gen.o SpMV_inspect.o

default: SpMV

SpMV: $(OBJS) 
	$(CC) -o $@ $^ $(LIB_PATHS) -lcudart 

%.o: %.cu
	$(NVCC) $(OPTFLAGS) $(INCL_PATHS) -c $< -o $@

%.o: %.cpp
	$(CC) $(OPTFLAGS) $(INCL_PATHS) -c $< -o $@

clean: 
	rm SpMV
	rm -f *.o
	rm -f *.linkinfo

