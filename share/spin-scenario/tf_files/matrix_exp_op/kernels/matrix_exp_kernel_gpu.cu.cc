// 2018, Cecilia Curreli

#if GOOGLE_CUDA

#define EIGEN_USE_GPU

#include "matrix_exp_op.h"
#include <cuda.h>
#include "tensorflow/core/util/cuda_kernel_helper.h"


//#include <device_launch_parameters.h>
//#include "third_party/eigen3/unsupported/Eigen/CXX11/Tensor"
//#include "/home/yan/spin-scenario/3rd-party/eigen/unsupported/Eigen/MatrixFunctions"
////#include "/home/yan/spin-scenario/3rd-party/eigen/Eigen"
//
//#include <cuda_runtime.h>
//#include <device_launch_parameters.h>
//#include <device_functions.h>
//#include <cuda_runtime_api.h>
////#include "external/cub_archive/cub/util_ptx.cuh"
//#include <cuComplex.h>
//#include <complex>
////#include <iostream.h>
#include <cstdlib> //cublas lib
#include <cublas_v2.h> //cublas lib


//Note: inside _global_ functions (functions on the GPU) we can't allocate new memory
// because tensorflow wants to monopolize the management of memory.
//This means, making matrix-matrix multiplication is complicated. We can't do A=A*B because we will overwrite A during the operation.
//If we want to have temporary tensor we can allocate them in the Functors.


namespace tensorflow {
namespace {

using CudaLaunchConfig = ::tensorflow::CudaLaunchConfig;

template <typename T>
__global__ void forward(CudaLaunchConfig cfg, T* __restrict__ Z, const int N,
                        const T* __restrict__ X, const T* __restrict__ Y,
                        const T bias) {

    //consider that the input matrices were probably Rowmajor before being flattened. Now they are arrays.
    //I suggest to transpose the input matrices before doing some calculations. maybe cublas can do that for us.

//H_n = H
//help = H_n
//fact=1
//int taylor_term = 6;
//for i=1 i<=taylor_term
//  fact= fact*i //this is safe since it's scalar mult
//  exp=exp+H_n/fact //this addition can be done in place. I suppose every matrix is initialized to zero automatically. Here for i=1 exp=0matrix
//  help=H_n
//  H_n = H*help
//

  // for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x
  // * gridDim.x) {
  for (int i : CudaGridRangeX(cfg.virtual_thread_count)) {
    Z[i] = X[i] + Y[i] + (T)bias;
  }
}

template <typename T>
__global__ void backward(CudaLaunchConfig cfg, const T* __restrict__ top_diff,
                         const int N, T* __restrict__ grad_matrixA,
                         T* __restrict__ grad_matrixB) {
  // for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += blockDim.x
  // * gridDim.x) {
  for (int i : CudaGridRangeX(cfg.virtual_thread_count)) {
    grad_matrixA[i] = top_diff[i];
    grad_matrixB[i] = top_diff[i];
  }
}

}  // anonymous namespace

namespace functor {

template <typename Dtype>
struct MatrixExpFunctor<GPUDevice, Dtype> {
  static void launch(::tensorflow::OpKernelContext* context, const Tensor& mA_,
                      Tensor* mOut_, float delta_t_) {

    const GPUDevice& d = context->eigen_gpu_device();
    //find optimal configuration for launching the Cuda Kernel. N,d are not tested.
    const int N = mA_.NumElements();
    ::tensorflow::CudaLaunchConfig cfg = ::tensorflow::GetCudaLaunchConfig(N, d);

    forward<Dtype><<<cfg.block_count, cfg.thread_per_block, 0, d.stream()>>>(
        cfg, mC_->flat<Dtype>().data(), mA_.NumElements(),
        mA_.flat<Dtype>().data(), mB_.flat<Dtype>().data(), bias);

    if (!d.ok()) {
      context->SetStatus(
          tensorflow::errors::Internal("Failed launching MatrixExp on GPU"));
    }
  }
};

template struct MatrixExpFunctor<GPUDevice, std::complex<float>>;
template struct MatrixExpFunctor<GPUDevice, std::complex<double>>;

template <typename Dtype>
struct MatrixExpGrad<GPUDevice, Dtype> {
  static void launch(::tensorflow::OpKernelContext* context,
                     const Tensor& topdiff_, Tensor* mA_,
                     Tensor* grad_mA_, float delta_t_) {
    const int N = topdiff_.NumElements();
    const GPUDevice& d = context->eigen_gpu_device();

    ::tensorflow::CudaLaunchConfig cfg =
        ::tensorflow::GetCudaLaunchConfig(N, d);

    // // optional reset gradients before running a kernel
    // cudaMemset(grad_mA_->flat<Dtype>().data(), 0, N * sizeof(Dtype));
    // cudaMemset(grad_mB_->flat<Dtype>().data(), 0, N * sizeof(Dtype));

    // backward<Dtype>
    // <<< cfg.block_count, cfg.thread_per_block, 0,
    // context->eigen_gpu_device().stream() >>> (
    //   cfg,
    //   topdiff_.flat<Dtype>().data(),
    //   topdiff_.NumElements(),
    //   grad_mA_->flat<Dtype>().data(),
    //   grad_mB_->flat<Dtype>().data());

    // this function could be useful in the future
    //cudaMemcpy(grad_mA_->flat<Dtype>().data(), topdiff_.flat<Dtype>().data(), N * sizeof(Dtype), cudaMemcpyDeviceToDevice);


    if (!d.ok()) {
      context->SetStatus(tensorflow::errors::Internal(
          "Failed launching MatrixExpGrad on GPU"));
    }
  }
};

template struct MatrixExpGrad<GPUDevice, std::complex<float>>;
template struct MatrixExpGrad<GPUDevice, std::complex<double>>;

}  // namespace functor
}  // namespace tensorflow

#endif  // GOOGLE_CUDA
