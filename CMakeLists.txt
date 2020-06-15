cmake_minimum_required(VERSION 2.8)
cmake_minimum_required(VERSION 2.8)
project(Spin-Scenario)
set(ssl_major 1)
set(ssl_minor 0)
set(ssl_patch 0)
set(ssl_version ${ssl_major}.${ssl_minor}.${ssl_patch})
message(STATUS "Configuring Spin-Scenario ${ssl_version}")
message(STATUS "GCC version: ${CMAKE_CXX_COMPILER_VERSION}") 

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if (WIN32)
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -Wall -MP -openmp -bigobj")
add_definitions(-DH5_BUILT_AS_DYNAMIC_LIB)
set(BOOST_ROOT "C:/local/boost_1_72_0")
set(Boost_USE_STATIC_LIBS ON) # only find static libs
#set(ArrayFire_DIR "C:/Program Files/ArrayFire/v3/cmake")
elseif (UNIX)
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall -fopenmp")
ENDIF ()

set(TENSORFLOW OFF)
#-------------path of 3rd party libraries-------------
# special libs.
find_package(Boost COMPONENTS system filesystem iostreams regex REQUIRED)
if (UNIX)
find_package(FFTW REQUIRED)
find_package(NLopt REQUIRED)
endif()
find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(PythonLibs REQUIRED)
#find_package(ArrayFire REQUIRED)
#find_package(CUDA)

set(EXTERN_LIB_ROOT ${PROJECT_SOURCE_DIR}/3rd-party)

# lib dirs.
#BUILD_SHARED_LIBS
if (WIN32)
set(FFTW_INCLUDES ${EXTERN_LIB_ROOT}/fftw)
set(FFTW_LIBRARIES "${EXTERN_LIB_ROOT}/fftw/libfftw3-3.lib")
set(NLOPT_INCLUDE_DIRS "${EXTERN_LIB_ROOT}/nlopt/build/src/api")
set(NLOPT_LIBRARIES "${EXTERN_LIB_ROOT}/nlopt/build/Release/nlopt.lib")
set(LUA_LIBRARIES "${EXTERN_LIB_ROOT}/lua/build/Release/lua53.lib") #5.3.4
set(LINENOISE_LIBRARIES "${EXTERN_LIB_ROOT}/linenoise-ng/build/Release/linenoise.lib")
set(YACAS_LIBRARIES "${EXTERN_LIB_ROOT}/yacas/cyacas/libyacas/build/Release/yacas.lib")
elseif (UNIX)
set(LUA_LIBRARIES "${EXTERN_LIB_ROOT}/lua/build/liblua53.so") #5.3.4
set(LINENOISE_LIBRARIES "${EXTERN_LIB_ROOT}/linenoise-ng/build/liblinenoise.so")
set(YACAS_LIBRARIES "${EXTERN_LIB_ROOT}/yacas/cyacas/libyacas/build/libyacas.so")
ENDIF ()

# tensorflow
if (UNIX)
if(TENSORFLOW)
set(TENSORFLOW_ROOT ${EXTERN_LIB_ROOT}/tensorflow)
set(TF_INCLUDE_DIRS "${TENSORFLOW_ROOT}" "${TENSORFLOW_ROOT}/bazel-genfiles" "${TENSORFLOW_ROOT}/bazel-tensorflow/external/protobuf_archive/src")
set(NEW_OP_ROOT ${PROJECT_SOURCE_DIR}/share/spin-scenario/tf_files)
set(TF_LIBRARIES ${TENSORFLOW_ROOT}/bazel-bin/tensorflow/libtensorflow_cc.so ${TENSORFLOW_ROOT}/bazel-bin/tensorflow/libtensorflow_framework.so)
endif()
endif()

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
        ${ArrayFire_INCLUDE_DIRS}
        ${TF_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${PYTHON_INCLUDE_DIRS}
        ${NEW_OP_ROOT}
        )


option(BUILD_SHARED_LIBS "build shared library" ON)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

#-------------ssl kernel lib-------------
set(_src_root "${PROJECT_SOURCE_DIR}/src")
file(
    GLOB_RECURSE _src_list 
    LIST_DIRECTORIES false
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${_src_root}/*.c*"
    "${_src_root}/*.h*")

add_library(ssl SHARED ${_src_list})

foreach(_source IN ITEMS ${_src_list})
    get_filename_component(_source_path "${_source}" PATH)
    string(REPLACE "/" "\\" _source_path_msvc "${_source_path}")
    source_group("${_source_path_msvc}" FILES "${_source}")
endforeach()
if(WIN32)
    set(SSL_LIBRARIES ${LUA_LIBRARIES} ${Boost_LIBRARIES} ${NLOPT_LIBRARIES} ${FFTW_LIBRARIES} ${LINENOISE_LIBRARIES} ${YACAS_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${PYTHON_LIBRARIES}  ${ArrayFire_LIBRARIES})
else()
    set(SSL_LIBRARIES ${TF_LIBRARIES} ${LUA_LIBRARIES} ${Boost_LIBRARIES} ${NLOPT_LIBRARIES} ${FFTW_LIBRARIES} ${LINENOISE_LIBRARIES} ${YACAS_LIBRARIES} ${HDF5_CXX_LIBRARIES} ${PYTHON_LIBRARIES})
endif()
if (WIN32)
target_link_libraries(ssl ${SSL_LIBRARIES})
elseif (UNIX)
target_link_libraries(ssl ${SSL_LIBRARIES} dl)
ENDIF (WIN32)
#set_target_properties(ssl PROPERTIES VERSION ${ssl_version} SOVERSION ${ssl_major})

add_executable(spin-scenario main.cpp)
target_link_libraries(spin-scenario ssl)

install(TARGETS spin-scenario ssl
RUNTIME DESTINATION bin
LIBRARY DESTINATION lib
CONFIGURATIONS Release)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/examples" DESTINATION .)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/config" DESTINATION share/spin-scenario)
install(FILES ${TF_LIBRARIES} ${LUA_LIBRARIES} ${LINENOISE_LIBRARIES} ${YACAS_LIBRARIES} DESTINATION lib)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE" "${CMAKE_CURRENT_SOURCE_DIR}/README.md" DESTINATION share/spin-scenario/doc)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/tf_files/matrix_exp_op/__init__.py" DESTINATION share/spin-scenario/tf_files/matrix_exp_op)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/tf_files/matrix_exp_op/matrix_exp_op.so" DESTINATION share/spin-scenario/tf_files/matrix_exp_op)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/share/spin-scenario/tf_files/graph_generator.py" DESTINATION share/spin-scenario/tf_files)

if (WIN32)
# for windows dll install
install(FILES "${EXTERN_LIB_ROOT}/fftw/libfftw3-3.dll" DESTINATION bin)
install(FILES "${EXTERN_LIB_ROOT}/linenoise-ng/build/Release/linenoise.dll" DESTINATION bin)
install(FILES "${EXTERN_LIB_ROOT}/lua/build/Release/lua53.dll" DESTINATION bin)
install(FILES "${EXTERN_LIB_ROOT}/nlopt/build/Release/nlopt.dll" DESTINATION bin)
install(FILES "${EXTERN_LIB_ROOT}/yacas/cyacas/libyacas/build/Release/yacas.dll" DESTINATION bin)
install(FILES "${HDF5_INCLUDE_DIRS}/../bin/hdf5.dll" DESTINATION bin)
install(FILES "${HDF5_INCLUDE_DIRS}/../bin/hdf5_cpp.dll" DESTINATION bin)
endif()

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
set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/docs/source/media/logo.ico")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}")
SET(CPACK_NSIS_MODIFY_PATH ON)
install (TARGETS spin-scenario LIBRARY DESTINATION lib ARCHIVE DESTINATION lib RUNTIME DESTINATION bin COMPONENT app)
if (UNIX)
set(CPACK_DEBIAN_PACKAGE_NAME, "Spin-Scenario")
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION " Spin-Scenario library, a flexible scripting environment for realistic MR simulations")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Yan Chang")
set(CPACK_DEBIAN_PACKAGE_VERSION "${ssl_version}")
set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libboost-filesystem-dev, libboost-system-dev, libboost-iostreams-dev, libboost-regex-dev, libfftw3-dev, libnlopt-dev, libnlopt-cxx-dev, libhdf5-dev, libpython3-dev, gnuplot")
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
endif()
include (CPack)
#cpack -C Release
