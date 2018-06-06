PROJECT(IMT_IRGN)

FIND_PACKAGE(CUDA)
if(NOT CUDA_FOUND)
 message("CUDA not found!")
endif(NOT CUDA_FOUND)

include_directories (${CUDA_INCLUDE_DIRS})
link_libraries (${CUDA_LIBRARIES})
link_libraries (${CUDA_CUBLAS_LIBRARIES})

# dependencies
FIND_PACKAGE(Qt4 REQUIRED)

# source definitions
INCLUDE_DIRECTORIES( ${CMAKE_SOURCE_DIR}
                     "${IMT_IRGN_SOURCE_DIR}/gui"
                     "${IMT_IRGN_BINARY_DIR}/gui"
                     "${IMT_IRGN_SOURCE_DIR}"
                     "${IMT_IRGN_BINARY_DIR}")


INCLUDE_DIRECTORIES( "${AGILE_SOURCE_DIR}/include"
                     "${AGILE_BINARY_DIR}/include"
                    )

###########################################################################
# set Agile Debug Mode
###########################################################################
ADD_DEFINITIONS(-DAGILE_DEBUG)

link_libraries(agile)

# add Qt modules here, the include will setup QT_LIBRARIES
SET(QT_USE_QTNETWORK TRUE)
SET(QT_USE_QTWEBKIT TRUE)

INCLUDE(${QT_USE_FILE})

# .cxx sources
SET(IMT_IRGN_SRCS_CXX

  gui/imt_mrgui.cpp
  gui/widgets/pathsetting.cpp
  gui/widgets/specialsetting.cpp
  gui/widgets/irgn_para.cpp
  gui/widgets/cycleread.cpp
  gui/widgets/fileblock.cpp

  #irgn/irgn_all.hpp
  imt_mrcmd.cpp
  main.cpp
)

# files which need to be moc'd by Qt
SET(IMT_IRGN_MOC_SRCS

  gui/imt_mrgui.h
  gui/widgets/pathsetting.h
  gui/widgets/specialsetting.h
  gui/widgets/irgn_para.h
  gui/widgets/cycleread.h
  #gui/widgets/fileblock.h
)

SET(IMT_IRGN_UI
  
  gui/imt_mrgui.ui
)

SET(QtApp_RCCS images.qrc)

# build ui_XXX files from the XML-style .ui files
QT4_WRAP_UI(IMT_IRGN_SRCS_CXX ${IMT_IRGN_UI})
# this moc's the above variable and appends to the cxx sources
QT4_WRAP_CPP(IMT_IRGN_SRCS_CXX ${IMT_IRGN_MOC_SRCS})
# generate rules for building source files from the resources
QT4_ADD_RESOURCES(QtApp_RCC_SRCS ${QtApp_RCCS})

ADD_EXECUTABLE(IMT_IRGN ${IMT_IRGN_SRCS_CXX} ${QtApp_RCC_SRCS})

TARGET_LINK_LIBRARIES(IMT_IRGN ${QT_LIBRARIES})

cuda_add_cufft_to_target (IMT_IRGN)

#copy executable to defined path
add_custom_command(TARGET IMT_IRGN 
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:IMT_IRGN> "${AGILE_SOURCE_DIR}/bin")


