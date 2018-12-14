# $Id: Makefile 13 2009-02-19 14:59:20Z freiberger $

# name of the final executable (without extension)
EXECUTABLE := dense_multiplication

# objects to link for the executable
OBJECTS := dense_multiplication.o

include ../../Make.common

# all: build the executable
all: $(EXECUTABLE)$(EXEEXT)

# run: build the executable and run it
run: $(EXECUTABLE)$(EXEEXT)
	@echo "[EXECUTE] "$<
	@./$(EXECUTABLE)$(EXEEXT)

# clean up
clean:
	-@rm -f $(OBJECTS) $(EXECUTABLE)$(EXEEXT)

# some debug output
mkdebug:
	@echo "CXX = "$(CXX)
	@echo "CXXFLAGS = "$(CXXFLAGS)
	@echo "NVCXX = "$(NVCXX)
	@echo "NVCXXFLAGS = "$(NVCXXFLAGS)
	@echo "LDFLAGS = "$(LDFLAGS)
	@echo "LIBS = "$(LIBS)
	@echo "OBJECTS = "$(OBJECTS)
	@echo 

# End of $Id: Makefile 13 2009-02-19 14:59:20Z freiberger $
