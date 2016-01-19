# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

###########################################################################
# everything needed to create the documentation
###########################################################################
find_package (Doxygen)
if (DOXYGEN_FOUND)
  # test for programs needed to create maths in documentation
  find_package (LATEX)
  if (NOT LATEX_COMPILER)
    message (STATUS "LaTeX was not found.")
  endif()
  if (NOT MAKEINDEX_COMPILER)
    message (STATUS "makeindex was not found.")
  endif()
  if (NOT DVIPS_CONVERTER)
    message (STATUS "dvips was not found.")
  endif()

  configure_file ("${AGILE_SOURCE_DIR}/doc/Doxyfile.in"
                  "${AGILE_BINARY_DIR}/doc/Doxyfile")
  add_custom_target (doc ${DOXYGEN_EXECUTABLE} "${AGILE_BINARY_DIR}/doc/Doxyfile"
                     WORKING_DIRECTORY ${AGILE_SOURCE_DIR})
endif()

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.
