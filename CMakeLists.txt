cmake_minimum_required(VERSION 2.8)
cmake_minimum_required(VERSION 2.8)
project(spin-scenario)
set(ssl_major 1)
set(ssl_minor 0)
set(ssl_patch 0)
set(ssl_version ${ssl_major}.${ssl_minor}.${ssl_patch})
message(STATUS "configuring spin-scenario ${ssl_version}")
message(STATUS "GCC version: ${CMAKE_CXX_COMPILER_VERSION}") 

set(CMAKE_BUILD_TYPE "Release") # Debug Release
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -std=c++14 -O3 -Wall -fopenmp -Wno-sign-compare")
SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#-------------path of 3rd party libraries-------------
# special libs.
find_package(Boost COMPONENTS system filesystem iostreams regex REQUIRED)
find_package(FFTW REQUIRED)
find_package(NLopt REQUIRED)
find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(PythonLibs REQUIRED)
#find_package(CUDA)

set(EXTERN_LIB_ROOT ${PROJECT_SOURCE_DIR}/3rd-party)

# lib dirs.
set(LUA_LIBRARIES "${EXTERN_LIB_ROOT}/lua/build/liblua53.so") #5.3.4
set(LINENOISE_LIBRARIES "${EXTERN_LIB_ROOT}/linenoise-ng/build/liblinenoise.so")
set(YACAS_LIBRARIES "${EXTERN_LIB_ROOT}/yacas/cyacas/libyacas/build/libyacas.so")

# tensorflow
set(TENSORFLOW_ROOT ${EXTERN_LIB_ROOT}/tensorflow)
set(TF_INCLUDE_DIRS "${TENSORFLOW_ROOT}" "${TENSORFLOW_ROOT}/bazel-genfiles" "${TENSORFLOW_ROOT}/bazel-tensorflow/external/protobuf_archive/src")
set(NEW_OP_ROOT ${PROJECT_SOURCE_DIR}/share/spin-scenario/tf_files)
set(TF_LIBRARIES ${TENSORFLOW_ROOT}/bazel-bin/tensorflow/libtensorflow_cc.so ${TENSORFLOW_ROOT}/bazel-bin/tensorflow/libtensorflow_framework.so)

#-------------ssl headers-------------
include_directories(${PROJECT_SOURCE_DIR}/src
        ${EXTERN_LIB_ROOT}/eigen
        ${EXTERN_LIB_ROOT}/gnuplot-iostream
        ${EXTERN_LIB_ROOT}/
        ${EXTERN_LIB_ROOT}/linenoise-ng/include
        ${EXTERN_LIB_ROOT}/yacas/cyacas/libyacas/include
        ${EXTERN_LIB_ROOT}/lua/src
        ${NLOPT_INCLUDE_DIRS}
        ${FFTW_INCLUDES}
        ${TF_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
        ${NEW_OP_ROOT}
        )


option(BUILD_SHARED_LIBS "build shared library" ON)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#-------------ssl kernel lib-------------
file(GLOB_RECURSE _src_list
        LIST_DIRECTORIES false
        RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/*.h" "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_library(ssl SHARED ${_src_list})
if(CUDA_FOUND)
    set(SSL_LIBRARIES ${TF_LIBRARIES} ${LUA_LIBRARIES} ${Boost_LIBRARIES} ${NLOPT_LIBRARIES} ${FFTW_LIBRARIES} ${LINENOISE_LIBRARIES} ${YACAS_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${PYTHON_LIBRARIES}  ${CUDA_LIBRARIES})
else()
    set(SSL_LIBRARIES ${TF_LIBRARIES} ${LUA_LIBRARIES} ${Boost_LIBRARIES} ${NLOPT_LIBRARIES} ${FFTW_LIBRARIES} ${LINENOISE_LIBRARIES} ${YACAS_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${PYTHON_LIBRARIES})
endif()

target_link_libraries(ssl ${SSL_LIBRARIES} dl)
#set_target_properties(ssl PROPERTIES VERSION ${ssl_version} SOVERSION ${ssl_major})

add_executable(spin-scenario main.cpp)
target_link_libraries(spin-scenario ssl)

install(TARGETS spin-scenario ssl
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
CONFIGURATIONS Release)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/config" DESTINATION share/spin-scenario)
install(FILES ${TF_LIBRARIES} ${LUA_LIBRARIES} ${LINENOISE_LIBRARIES} ${YACAS_LIBRARIES} DESTINATION lib)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE" "${CMAKE_CURRENT_SOURCE_DIR}/README.md" DESTINATION share/spin-scenario/doc)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/tf_files/matrix_exp_op/__init__.py" DESTINATION share/spin-scenario/tf_files/matrix_exp_op)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/tf_files/matrix_exp_op/matrix_exp_op.so" DESTINATION share/spin-scenario/tf_files/matrix_exp_op)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/tf_files/graph_generator.py" DESTINATION share/spin-scenario/tf_files)

include(InstallRequiredSystemLibraries)
# Create Ubuntu package
set(CPACK_PACKAGE_VENDOR "The Spin-Scenario Team")
set(CPACK_PACKAGE_CONTACT "changy@sibet.ac.cn")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${ssl_major}")
set(CPACK_PACKAGE_VERSION_MINOR "${ssl_minor}")
set(CPACK_PACKAGE_VERSION_PATCH "${ssl_patch}")
set(CPACK_PACKAGE_RELEASE 1)

install (TARGETS spin-scenario LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin COMPONENT app)
set(CPACK_DEBIAN_PACKAGE_NAME, "Spin-Scenario")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION " Spin-Scenario library, a flexible scripting environment for realistic MR simulations")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Yan Chang")
set(CPACK_DEBIAN_PACKAGE_VERSION "${ssl_version}")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-filesystem-dev, libboost-system-dev, libboost-iostreams-dev, libboost-regex-dev, libfftw3-dev, libnlopt-dev, libhdf5-dev, libpython3-dev, gnuplot")
set(CPACK_GENERATOR "DEB")

# RPM package
#set(CPACK_RPM_PACKAGE_NAME "Spin-Scenario")
#set(CPACK_RPM_PACKAGE_GROUP "The Spin-Scenario Team")
#set(CPACK_RPM_PACKAGE_URL "https://github.com/spin-scenario")
#set(CPACK_RPM_PACKAGE_VERSION "${ssl_version}")
#set(CPACK_RPM_PACKAGE_RELEASE 1)
#set(CPACK_RPM_PACKAGE_VENDOR "The Spin-Scenario Team")
#set(CPACK_RPM_PACKAGE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
#set(CPACK_RPM_PACKAGE_DESCRIPTION " Spin-Scenario library, a flexible scripting environment for realistic MR simulations")
#set(CPACK_RPM_PACKAGE_REQUIRES "boost-devel, fftw-devel, hdf5-devel, NLopt-devel, gnuplot") #python3
#set(CPACK_GENERATOR "RPM")

include (CPack)
#cpack -C Release
