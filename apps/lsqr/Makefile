# $Id: Make.test 128 2009-03-24 13:41:21Z freiberger $

# name of the final executable (without extension)
EXECUTABLE := lsqr

# objects to link for the executable
SOURCES := lsqr.cpp

include ../../Make.common

# all: build the executable
all: $(EXECUTABLE)$(MODEEXEEXT)

# run: build the executable and run it
run: $(EXECUTABLE)$(MODEEXEEXT)
	@clear
	@echo "[EXECUTE] "$<
	@./$(EXECUTABLE)$(MODEEXEEXT) ./data/csrmatrix_A_1048576x1048576.bin ./data/csrmatrix_AT_1048576x1048576.bin ./data/vector_kspace_data_1048576\(1024x1024\).bin ./x_result.bin

# clean up
clean:
	-@rm -f $(patsubst %.cu,%$(CUOBJEXT),$(patsubst %.cpp,%$(OBJEXT),$(SOURCES)))
	-@rm -f $(patsubst %.cu,%-dbg$(CUOBJEXT),$(patsubst %.cpp,%-dbg$(OBJEXT),$(SOURCES)))
	-@rm -f $(EXECUTABLE)$(EXEEXT) $(EXECUTABLE)-dbg$(EXEEXT)

# some debug output
make_debug:
	@echo "CXX = "$(CXX)
	@echo "CXXFLAGS = "$(CXXFLAGS)
	@echo "NVCXX = "$(NVCXX)
	@echo "NVCXXFLAGS = "$(NVCXXFLAGS)
	@echo "LDFLAGS = "$(LDFLAGS)
	@echo "LIBS = "$(LIBS)
	@echo "SOURCES = "$(SOURCES)
	@echo 

# End of $Id: Make.test 128 2009-03-24 13:41:21Z freiberger $.
