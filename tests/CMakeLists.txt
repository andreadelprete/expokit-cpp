SET(DISABLE_TESTS OFF)

MACRO(ADD_TEST_CFLAGS target flag)
  SET_PROPERTY(TARGET ${target} APPEND_STRING PROPERTY COMPILE_FLAGS " ${flag}")
ENDMACRO(ADD_TEST_CFLAGS)

MACRO(ADD_TESTCASE NAME PKGS)
  ADD_UNIT_TEST(${NAME} ${NAME})

  SET_TARGET_PROPERTIES(${NAME} PROPERTIES LINKER_LANGUAGE CXX)
  #ADD_TEST_CFLAGS(${NAME} "-DBOOST_TEST_DYN_LINK") 
  SET(MODULE_NAME "${NAME}Test")
  STRING(REGEX REPLACE "-" "_" MODULE_NAME ${MODULE_NAME})
  ADD_TEST_CFLAGS(${NAME} "-DBOOST_TEST_MODULE=${MODULE_NAME}") 

  FOREACH(PKG ${PKGS})
    PKG_CONFIG_USE_DEPENDENCY(${NAME} ${PKG})
  ENDFOREACH(PKG)

  TARGET_LINK_LIBRARIES(${NAME} ${PROJECT_NAME})
  #TARGET_LINK_LIBRARIES(${NAME} ${LAPACK_LIBRARIES}) # uncomment this line to use LAPACK
ENDMACRO(ADD_TESTCASE)

# --- RULES -------------------------------------------------------------------
#ADD_TESTCASE(test_dense_small "f2c")
#ADD_TESTCASE(test_dense_large "f2c")
#ADD_TESTCASE(test_small       "eigen3")

ADD_TESTCASE(testNoMalloc "eigen3")
ADD_TESTCASE(testLDS "eigen3")
ADD_TESTCASE(testMatrixExp "eigen3")
ADD_TESTCASE(testBalancing "eigen3")
ADD_TESTCASE(testBothIntegrals "eigen3")



#[[
]]



