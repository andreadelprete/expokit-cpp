
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Delta performance evaluation
add_executable(runDelta runDelta.cpp)
target_include_directories(runDelta PRIVATE ../include)
target_sources(runDelta PRIVATE ../src/stop-watch.cpp)
target_link_libraries (runDelta Eigen3::Eigen)

# Simple example of basic usage
add_executable(runExample runExample.cpp)
target_include_directories(runExample PRIVATE ../include)
target_sources(runExample PRIVATE ../src/stop-watch.cpp)
target_link_libraries (runExample Eigen3::Eigen)

# Example of TSV loader
add_executable(runTSVExample runTSVExample.cpp)
target_include_directories(runTSVExample PRIVATE ../include)
target_sources(runTSVExample PRIVATE ../src/stop-watch.cpp)
target_link_libraries (runTSVExample Eigen3::Eigen)

# Delta and TS performance evaluation
add_executable(runBench runBench.cpp)
target_include_directories(runBench PRIVATE ../include)
target_sources(runBench PRIVATE ../src/stop-watch.cpp)
target_link_libraries (runBench Eigen3::Eigen)

