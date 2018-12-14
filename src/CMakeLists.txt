# $Id: CMakeLists.txt 493 2011-07-26 10:36:52Z huber $

###########################################################################
# This is a list of sources to be compile into the library.
###########################################################################
set (SOURCES
     ${CMAKE_CURRENT_SOURCE_DIR}/gpu_environment.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/quadtree.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/io/io_utils.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/io/readSiemensVD11.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/io/dicom.cpp
    )

###########################################################################
# This is a list of NVCC-sources to be compiled into the library.
###########################################################################
set (NVCC_SOURCES
     ${CMAKE_CURRENT_SOURCE_DIR}/gpu_communicator.cu
     ${CMAKE_CURRENT_SOURCE_DIR}/gpu_memory.cu
    )

# Function to determine if a type is complex
function (iscomplex ReturnVariable Type)
  string (REGEX MATCH "complex" MATCHINGSTRING ${Type})
  string (COMPARE EQUAL "complex" "${MATCHINGSTRING}" TypeIsComplex)
  set (${ReturnVariable} ${TypeIsComplex} PARENT_SCOPE)
endfunction (iscomplex)

# Macro to automatically generate files with specialized type (1 different type).
macro (generate1 Basename Extension Type1)
  # Set the types that will be replaced in the file.
  set (TType1 ${Type1})
  iscomplex (TType1IsComplex ${Type1})
  # Clean the input types so that we can use them to generate unique file names.
  string (REGEX REPLACE "[^a-zA-Z0-9_]" "" FType1 "${Type1}")
  # This is a warning message included in the automatically generated files.
  set (AutoWarning "/* This file was generated automatically by CMake. You have to modify '${CMAKE_CURRENT_SOURCE_DIR}/${Basename}.${Extension}.in' if you want to make changes. */")
  configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/${Basename}.${Extension}.in"
                  "${CMAKE_CURRENT_BINARY_DIR}/auto/${Basename}_${FType1}.${Extension}")
  # add the file to the list of files to be compiled
  list (APPEND AUTO_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/auto/${Basename}_${FType1}.${Extension}")
  unset (AutoWarning)
  unset (TType1)
  unset (TType1IsComplex)
endmacro (generate1)

# Macro to automatically generate files with specialized types (2 different types).
macro (generate2 Basename Extension Type1 Type2)
  # Set the types that will be replaced in the file.
  set (TType1 ${Type1})
  set (TType2 ${Type2})
  iscomplex (TType1IsComplex ${Type1})
  iscomplex (TType2IsComplex ${Type2})
  # Clean the input types so that we can use them to generate unique file names.
  string (REGEX REPLACE "[^a-zA-Z0-9_]" "" FType1 "${Type1}")
  string (REGEX REPLACE "[^a-zA-Z0-9_]" "" FType2 "${Type2}")
  # This is a warning message included in the automatically generated files.
  set (AutoWarning "/* This file was generated automatically by CMake. You have to modify '${CMAKE_CURRENT_SOURCE_DIR}/${Basename}.${Extension}.in' if you want to make changes. */")
  configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/${Basename}.${Extension}.in"
                  "${CMAKE_CURRENT_BINARY_DIR}/auto/${Basename}_${FType1}_${FType2}.${Extension}")
  # add the file to the list of files to be compiled
  list (APPEND AUTO_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/auto/${Basename}_${FType1}_${FType2}.${Extension}")
  unset (AutoWarning)
  unset (TType1)
  unset (TType2)
  unset (TType1IsComplex)
  unset (TType2IsComplex)
endmacro (generate2)

# Macro to automatically generate files with specialized types (3 different types).
macro (generate3 Basename Extension Type1 Type2 Type3)
  # Set the types that will be replaced in the file.
  set (TType1 ${Type1})
  set (TType2 ${Type2})
  set (TType3 ${Type3})
  iscomplex (TType1IsComplex ${Type1})
  iscomplex (TType2IsComplex ${Type2})
  iscomplex (TType3IsComplex ${Type3})
  string (COMPARE EQUAL "${Type2}" "${Type3}" TType2EqualTType3)
  # Clean the input types so that we can use them to generate unique file names.
  string (REGEX REPLACE "[^a-zA-Z0-9_]" "" FType1 "${Type1}")
  string (REGEX REPLACE "[^a-zA-Z0-9_]" "" FType2 "${Type2}")
  string (REGEX REPLACE "[^a-zA-Z0-9_]" "" FType3 "${Type3}")
  # This is a warning message included in the automatically generated files.
  set (AutoWarning "/* This file was generated automatically by CMake. You have to modify '${CMAKE_CURRENT_SOURCE_DIR}/${Basename}.${Extension}.in' if you want to make changes. */")
  configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/${Basename}.${Extension}.in"
                  "${CMAKE_CURRENT_BINARY_DIR}/auto/${Basename}_${FType1}_${FType2}_${FType3}.${Extension}")
  # add the file to the list of files to be compiled
  list (APPEND AUTO_SOURCES "${CMAKE_CURRENT_BINARY_DIR}/auto/${Basename}_${FType1}_${FType2}_${FType3}.${Extension}")
  unset (AutoWarning)
  unset (TType1)
  unset (TType2)
  unset (TType3)
  unset (TType1IsComplex)
  unset (TType2IsComplex)
  unset (TType3IsComplex)
  unset (TType2EqualTType3)
endmacro (generate3)

###########################################################################
# Generate all the files in the "auto"-subdirectory
###########################################################################

generate2 (gpu_cs_matrix cu float               float)
generate2 (gpu_cs_matrix cu std::complex<float> float)
generate2 (gpu_cs_matrix cu float               std::complex<float>)
generate2 (gpu_cs_matrix cu std::complex<float> std::complex<float>)

generate2 (gpu_matrix_pitched cu float               float)
generate2 (gpu_matrix_pitched cu std::complex<float> float)
generate2 (gpu_matrix_pitched cu float               std::complex<float>)
generate2 (gpu_matrix_pitched cu std::complex<float> std::complex<float>)

generate1 (gpu_matrix cpp float)
generate1 (gpu_matrix cpp std::complex<float>)


# integer types
generate3 (gpu_vector cu "unsigned char"     "unsigned char"     "unsigned char")
generate3 (gpu_vector cu int                 int                 int)
generate3 (gpu_vector cu unsigned            unsigned            unsigned)
# floating point types
generate3 (gpu_vector cu float               float               float)
generate3 (gpu_vector cu std::complex<float> float               float)
generate3 (gpu_vector cu float               std::complex<float> float)
generate3 (gpu_vector cu std::complex<float> std::complex<float> float)
generate3 (gpu_vector cu float               float               std::complex<float>)
generate3 (gpu_vector cu std::complex<float> float               std::complex<float>)
generate3 (gpu_vector cu float               std::complex<float> std::complex<float>)
generate3 (gpu_vector cu std::complex<float> std::complex<float> std::complex<float>)

generate1 (calc/fft cpp std::complex<float>)
generate1 (calc/genkspacefov cpp std::complex<float>)
generate2 (calc/irgn cpp std::complex<float> float)
generate2 (calc/l2solve cpp std::complex<float> float)
generate2 (calc/tvsolve cpp std::complex<float> float)
generate2 (calc/tgvsolve cpp std::complex<float> float)
generate2 (calc/postprocess cpp std::complex<float> float)

# Double types, only supported when Cuda Architecture 200 is used
if(ENABLE_GPU_DOUBLE)
    generate1 (gpu_matrix cpp double)
    generate1 (gpu_matrix cpp std::complex<double>)

    generate2 (gpu_matrix_pitched cu double              double)
    generate2 (gpu_matrix_pitched cu std::complex<double> double)
    generate2 (gpu_matrix_pitched cu double               std::complex<double>)
    generate2 (gpu_matrix_pitched cu std::complex<double> std::complex<double>)

    generate3 (gpu_vector cu double               double               double)
    generate3 (gpu_vector cu std::complex<double> double               double)
    generate3 (gpu_vector cu double               std::complex<double> double)
    generate3 (gpu_vector cu std::complex<double> std::complex<double> double)
    generate3 (gpu_vector cu double               double               std::complex<double>)
    generate3 (gpu_vector cu std::complex<double> double               std::complex<double>)
    generate3 (gpu_vector cu double               std::complex<double> std::complex<double>)
    generate3 (gpu_vector cu std::complex<double> std::complex<double> std::complex<double>)

    generate1 (calc/fft cpp std::complex<double>)
    generate1 (calc/genkspacefov cpp std::complex<double>)
    generate2 (calc/irgn cpp std::complex<double> double)
    generate2 (calc/l2solve cpp std::complex<double> double)
    generate2 (calc/tvsolve cpp std::complex<double> double)
    generate2 (calc/tgvsolve cpp std::complex<double> double)
    generate2 (calc/postprocess cpp std::complex<double> double)
endif(ENABLE_GPU_DOUBLE)

# Add include directory for the header files. Note that we also have to include
# the binary directory because the CONFIGURE_FILE command creates the gpu_config.hpp
# in there. Also the current source directory is needed for the *.ipp files.
include_directories ("${AGILE_SOURCE_DIR}/include"
                     "${AGILE_BINARY_DIR}/include"
                     "${CMAKE_CURRENT_SOURCE_DIR}"
                    )

cuda_include_directories ("${AGILE_SOURCE_DIR}/include"
                          "${AGILE_BINARY_DIR}/include"
                          "${CMAKE_CURRENT_SOURCE_DIR}"
                         )

# Create a library
cuda_add_library (agile
                  ${NVCC_SOURCES} ${AUTO_SOURCES} ${SOURCES}
                  STATIC
                  # OPTIONS "-DSTUFF=blah blah"
                  # RELEASE --use_fast_math -DNDEBUG
                  # DEBUG -g -DDEBUG
                 )

# Install the library
install (TARGETS agile
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib
        )

# End of $Id: CMakeLists.txt 493 2011-07-26 10:36:52Z huber $.
