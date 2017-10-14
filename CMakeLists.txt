cmake_minimum_required(VERSION 3.3.0)
project(PPA_Assembler)

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/release/)
#set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib/)

set(CMAKE_C_COMPILER   "/data/opt/brew/bin/mpicc")
set(CMAKE_CXX_COMPILER "/data/opt/brew/bin/mpic++")
set(CMAKE_BUILD_TYPE Release CACHE STRING "set build type to release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
set(COMMON_LINK_LIBS "-lhdfs -lpthread") #-ljvm 

### mpicc&mpic++ will fix MPI_INCLUDE&MPI_LIBRARY themselves
#set(MPI_INCLUDE_DIR "/data/opt/brew/include")
#set(MPI_LIBRARY_DIR "/data/opt/brew/lib")

set(HDFS_INCLUDE_DIR "/data/opt/hadoop-2.6.0/include/")
set(HDFS_LIBRARY_DIR "/data/opt/hadoop-2.6.0/lib/native/")

set(PPA_Assembler_EXTERNAL_INCLUDES  ${HDFS_INCLUDE_DIR})
set(PPA_Assembler_EXTERNAL_LIBRARIES ${HDFS_LIBRARY_DIR})


add_subdirectory(utils)

add_subdirectory(put)
add_subdirectory(example)
