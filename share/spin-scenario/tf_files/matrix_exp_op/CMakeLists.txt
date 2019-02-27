#2018, Cecilia Curreli
cmake_minimum_required( VERSION 2.8 )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}../../../../cmake)

project(TFCustomOp)

find_package(TensorFlow 10 REQUIRED) #tensorflow 9 is not fully tested
#find_package(CUDA 9 REQUIRED)

# set necessary flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SSE_FLAGS} -march=native -fopenmp -D_GLIBCXX_USE_CXX11_ABI=${TensorFlow_ABI}")
set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -fPIC --shared -D_GLIBCXX_USE_CXX11_ABI=${TensorFlow_ABI}" )
#pay attention to --gpu-architecture=sm_35 if you have a different gpu architecture
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -std=c++11  --expt-relaxed-constexpr -D GOOGLE_CUDA=1 --gpu-architecture=sm_35 -D_GLIBCXX_USE_CXX11_ABI=${TensorFlow_ABI}" )
#-gencode arch=compute_35
include_directories(SYSTEM "${CUDA_INCLUDE_DIRS}/../../")
include_directories(SYSTEM ${TensorFlow_INCLUDE_DIR})
include_directories(SYSTEM "kernels")

if(CUDA_FOUND)
	add_tensorflow_cpu_operation("matrix_exp")
#	add_tensorflow_gpu_operation("matrix_exp") #gpu operation is not ready yet
else()
 	add_tensorflow_cpu_operation("matrix_exp")
endif()



#arch=compute_${CMAKE_MATCH_2},code=sm_${CMAKE_MATCH_1})

