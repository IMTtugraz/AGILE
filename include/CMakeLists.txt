# $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $

ADD_SUBDIRECTORY (agile)

# Install the sub-directory "agile" into the "include"-directory relative to the
# CMAKE_INSTALL_PREFIX.
# have a look at http://www.cmake.org/cmake/help/cmake2.6docs.html#command:install
INSTALL (DIRECTORY agile
         DESTINATION include
         FILES_MATCHING
         PATTERN "*.hpp"
         PATTERN "*.h" 
         PATTERN ".svn" EXCLUDE)

# End of $Id: CMakeLists.txt 476 2011-06-16 08:54:14Z freiberger $.
