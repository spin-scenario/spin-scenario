// 2018, Cecilia Curreli

#ifndef MATRIX_EXP_KERNELS_MATRIX_EXP_OP_H_
#define MATRIX_EXP_KERNELS_MATRIX_EXP_OP_H_

#include "tensorflow/core/framework/op_kernel.h"

namespace tensorflow {

typedef Eigen::ThreadPoolDevice CPUDevice;
typedef Eigen::GpuDevice GPUDevice;

namespace functor {

template <typename Device, typename Dtype>
struct MatrixExpFunctor {
  static void launch(::tensorflow::OpKernelContext* ctx, const Tensor& X,
                      Tensor* Z, float delta_t);
};

template <typename Device, typename Dtype>
struct MatrixExpGrad {
  static void launch(::tensorflow::OpKernelContext* ctx,  const Tensor& topdiff_, const Tensor& X,
                     Tensor* grad_A, float delta_t);
};

}  // namespace functor
}  // namespace tensorflow

#endif  // MATRIX_EXP_KERNELS_MATRIX_EXP_OP_H_
