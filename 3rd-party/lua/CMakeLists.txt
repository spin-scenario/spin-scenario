cmake_minimum_required(VERSION 2.8)
project(lua53)
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -std=c++17")
file(GLOB_RECURSE _src_list
        LIST_DIRECTORIES false
        RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/*.h" "${PROJECT_SOURCE_DIR}/src/*.c")
add_library(lua53 SHARED ${_src_list})
target_link_libraries(lua53)
