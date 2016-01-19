# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

# List of source files which have to be compiled for the unit tests.
set (SOURCES gpu_vector_test.cpp
             #large_gpu_vector_test.cpp
             gpu_matrix_test.cpp
             #large_gpu_matrix_test.cpp
             gpu_cs_matrix_test.cpp
             test_main.cpp)

include_directories ("${AGILE_SOURCE_DIR}/include"
                     "${AGILE_BINARY_DIR}/include"
                    )

cuda_include_directories ("${AGILE_SOURCE_DIR}/include"
                          "${AGILE_BINARY_DIR}/include"
                         )

# We need the Boost library for their unit test framework.
find_package( Boost 1.36.0 COMPONENTS unit_test_framework )
if (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
  # The unittest executable.
  add_executable(unittest ${SOURCES})
  target_link_libraries(unittest ${Boost_LIBRARIES} agile)
  # A custom target to run the unittest.
  add_custom_target(unittest-run
                    unittest --report_format=XML --report_level=detailed 2> unittestreport.xml
                    DEPENDS unittest
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.

