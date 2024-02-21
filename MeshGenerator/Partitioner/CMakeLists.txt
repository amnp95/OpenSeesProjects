cmake_minimum_required(VERSION 3.0)
project(Partitioner)

# get GKLip path from user using options
option(GKLIB_PATH "Path to GKlib" "")
if(GKLIB_PATH)
    include_directories(${GKLIB_PATH}/include)
    link_directories(${GKLIB_PATH}/lib)
else()
    message(FATAL_ERROR "GKlib path not provided")
endif()

# get METIS path from user using options
option(METIS_PATH "Path to METIS" "")
if(METIS_PATH)
    include_directories(${METIS_PATH}/include)
    link_directories(${METIS_PATH}/lib)
else()
    message(FATAL_ERROR "METIS path not provided")
endif()

# shared library
add_library(Partitioner SHARED Partition.c)
target_link_libraries(Partitioner metis GKlib)


# install the library in the major folder
install(TARGETS Partitioner DESTINATION ${CMAKE_SOURCE_DIR}/../lib)



