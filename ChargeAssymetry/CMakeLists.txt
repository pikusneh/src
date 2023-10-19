
cmake_minimum_required(VERSION 3.0)

project(TopAnalysis)

#---Create a shared library with geneated dictionary
set (CMAKE_CXX_FLAGS_INIT "-g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

include_directories("${PROJECT_SOURCE_DIR};$ENV{ROOTSYS}/include/")

include_directories("${PROJECT_SOURCE_DIR}/include")

add_executable(chargeAssym  ElecMuon.cpp)

aux_source_directory(src top_sources)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS RIO Net)

add_library(chargeAssym SHARED ${top_sources})

target_link_libraries(chargeAssym ${ROOT_LIBRARIES})


