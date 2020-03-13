MACRO(ADD_RUN NAME PKGS)
  add_executable(${NAME} ${NAME}.cpp)

  SET_TARGET_PROPERTIES(${NAME} PROPERTIES LINKER_LANGUAGE CXX)
  SET(MODULE_NAME "${NAME}Test")
  STRING(REGEX REPLACE "-" "_" MODULE_NAME ${MODULE_NAME})

  FOREACH(PKG ${PKGS})
    PKG_CONFIG_USE_DEPENDENCY(${NAME} ${PKG})
  ENDFOREACH(PKG)

  TARGET_LINK_LIBRARIES(${NAME} ${PROJECT_NAME})
ENDMACRO(ADD_RUN)

ADD_RUN(runDelta "eigen3")
ADD_RUN(runBench "eigen3")
ADD_RUN(runExample "eigen3")
ADD_RUN(runTSVExample "eigen3")
ADD_RUN(runFindMinSquarings "eigen3")
ADD_RUN(runFindVecSquarings "eigen3")


ADD_RUN(runBalancing "eigen3")






