# List of source files which have to be compiled for this test.

#SET(GDCM_DIR ${AGILE_SOURCE_DIR}/../gdcm-2.2.0/gdcm-build)
#FIND_PACKAGE(GDCM REQUIRED)
#INCLUDE(${GDCM_USE_FILE})

#find_package(GDCM)

#set (SOURCES irgnmain.cpp)
set (SOURCES cufft_test.cpp)



#SET(GDCM_LIBRARIES gdcmDICT gdcmDSED gdcmIOD
#    gdcmjpeg8 gdcmjpeg12 gdcmjpeg16 gdcmMSFF gdcmMEXD
#    CACHE STRING "GDCM libraries")

#MESSAGE(STATUS "Cmake GDCM: using includes at: ${GDCM_INCLUDE_DIRS}")
#MESSAGE(STATUS "Cmake GDCM: using libraries at: ${GDCM_LIBRARY_DIRS}")
#MESSAGE(STATUS "Cmake GDCM: using libraries: ${GDCM_LIBRARYRIES}")

#TARGET_LINK_LIBRARIES(irgnmain ${GDCM_LIBRARIES})

# The resulting executable.
add_executable(irgnmain ${SOURCES})

#target_link_libraries(irgnmain gdcmMSFF gdcmCommon gdcmDICT gdcmDSED gdcmIOD)
#target_link_libraries(irgnmain gdcmCommon gdcmDSED)

cuda_add_cufft_to_target (irgnmain)












# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

#include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List of source files which have to be compiled for this test.
#set (SOURCES1 irgnmain.cpp irgn.cpp l2solve.hpp tvsolve.hpp tgvsolve.hpp)

# The resulting executable.
#cuda_add_executable(IRGN ${SOURCES1})
#cuda_add_cufft_to_target (IRGN)


# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.
