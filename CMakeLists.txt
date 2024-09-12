cmake_minimum_required(VERSION 3.1)
project(WTCT 
  VERSION 0.0.1
  DESCRIPTION "HLL: HyperLogLog implementation"
  LANGUAGES CXX
)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(
  ${CMAKE_HOME_DIRECTORY}/include
  ${CMAKE_HOME_DIRECTORY}/external/smhasher/src
  ${STXXL_INCLUDE_DIRS}
)

add_subdirectory(external/smhasher/src)

add_executable(main src/main.cpp)
target_link_libraries(main SMHasherSupport)



