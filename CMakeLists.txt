CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
#INCLUDE(cmake/boost.cmake)
#INCLUDE(cmake/eigen.cmake)
#INCLUDE(cmake/python.cmake)
INCLUDE(cmake/lapack.cmake)
INCLUDE(cmake/ide.cmake)
INCLUDE(cmake/apple.cmake)
include(cmake/test.cmake)

SET(PROJECT_NAME expokit)
SET(PROJECT_DESCRIPTION "")
SET(PROJECT_URL "")
OPTION(INSTALL_DOCUMENTATION "Generate and install the documentation" OFF)

SET(DOXYGEN_USE_MATHJAX YES)

# Handle APPLE Cmake policy
IF(APPLE)
  APPLY_DEFAULT_APPLE_CONFIGURATION()
ENDIF(APPLE)

# Disable -Werror on Unix for now.
SET(CXX_DISABLE_WERROR True)
SET(CMAKE_VERBOSE_MAKEFILE True)

SETUP_PROJECT()

IF(WIN32)
  SET(LINK copy_if_different)
ELSE(WIN32)
  SET(LINK create_symlink)
ENDIF(WIN32)


# ----------------------------------------------------
# --- DEPENDENCIES -----------------------------------
# ----------------------------------------------------

ADD_DEFINITIONS(-DEIGEN_RUNTIME_NO_MALLOC)
ADD_DEFINITIONS(-DPROFILERYES)

# uncomment these lines to use LAPACK
#SEARCH_FOR_LAPACK()
#INCLUDE_DIRECTORIES(SYSTEM ${LAPACK_INCLUDE_DIRS})

ADD_REQUIRED_DEPENDENCY("f2c")
ADD_REQUIRED_DEPENDENCY("eigen3")

# Path to boost headers
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# ----------------------------------------------------
# --- INCLUDE ----------------------------------------
# ----------------------------------------------------
SET(HEADERS 
    include/clock.h
    include/dgchbv.h
    include/dgexpv.h
    include/dgpadm.h
    include/expokit.h
    include/MatrixExponential.hpp
    include/LDSUtility.hpp
    include/LDS2OrderUtility.hpp
    include/utils/stop-watch.h
    include/utils/Stdafx.h
)

#LIST(REMOVE_DUPLICATES HEADERS)

SET(HEADERS_FULL_PATH "")
FOREACH(header ${HEADERS})
  LIST(APPEND HEADERS_FULL_PATH "${CMAKE_SOURCE_DIR}/${header}")
  GET_FILENAME_COMPONENT(headerName ${header} NAME)
  GET_FILENAME_COMPONENT(headerPath ${header} PATH)
  INSTALL(FILES ${${PROJECT_NAME}_SOURCE_DIR}/${header}
    DESTINATION ${CMAKE_INSTALL_PREFIX}/${headerPath}
          PERMISSIONS OWNER_READ GROUP_READ WORLD_READ OWNER_WRITE)
ENDFOREACH(header ${HEADERS})

# --- MAIN LIBRARY -------------------------------------------------------------
ADD_SUBDIRECTORY(src)

# --- TESTS -------------------------------------------------------------
ADD_SUBDIRECTORY(tests)


# --- PACKAGING ----------------------------------------------------------------
PKG_CONFIG_APPEND_LIBS(${PROJECT_NAME})

SETUP_PROJECT_FINALIZE()

add_subdirectory(run) 

