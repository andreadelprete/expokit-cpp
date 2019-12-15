# --- MAIN LIBRARY -------------------------------------------------------------
SET(LIBRARY_NAME ${PROJECT_NAME})

SET(${LIBRARY_NAME}_SOURCES
    blas.c   
    dgchbv.c  
    dgexpv.c  
    dgpadm.c  
#    dmexpv.c  
#    dschbv.c  
#    dspadm.c
#    loadhb.c  
#    zgchbv.c  
#    zgexpv.c  
#    zgpadm.c  
#    zhexpv.c  
#    zhphiv.c
    clock.c  
#    dgcnvr.c  
#    dgmatv.c  
    dgphiv.c  
#    dnchbv.c  
#    dsexpv.c  
#    dsphiv.c  
    lapack.c   
    mataid.c
#    zgcnvr.c  
#    zgmatv.c  
#    zgphiv.c  
#    zhpadm.c  
#    znchbv.c
    expokit.cpp
  )

ADD_LIBRARY(${LIBRARY_NAME} SHARED ${HEADERS_FULL_PATH} ${${LIBRARY_NAME}_SOURCES})
SET_TARGET_PROPERTIES(${LIBRARY_NAME} PROPERTIES LINKER_LANGUAGE CXX)
ADD_HEADER_GROUP(HEADERS_FULL_PATH)
ADD_SOURCE_GROUP(${LIBRARY_NAME}_SOURCES)

PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} f2c)
PKG_CONFIG_USE_DEPENDENCY(${LIBRARY_NAME} eigen3)

INSTALL(TARGETS ${PROJECT_NAME} DESTINATION lib)

