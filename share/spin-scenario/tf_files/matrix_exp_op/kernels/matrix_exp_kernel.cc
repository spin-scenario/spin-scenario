// 2018, Cecilia Curreli

#include <cstring>

#include "matrix_exp_op.h"
#include "tensorflow/core/framework/op.h"

#include <unsupported/Eigen/MatrixFunctions>
#include <complex>


namespace tensorflow {
namespace functor {

template <typename Dtype>
struct MatrixExpFunctor<CPUDevice, Dtype> {
  static void launch(::tensorflow::OpKernelContext* ctx, const Tensor& Xt,
                      Tensor* Zt, float delta_t_) {

    auto X_dims = (Xt.shape()).dim_sizes();

      const auto mapIn = Eigen::Map<const Eigen::Matrix<
      Dtype,           /* scalar element type */
      Eigen::Dynamic,  /* num_rows is a run-time value */
      Eigen::Dynamic,  /* num_cols is a run-time value */
      Eigen::RowMajor  >>(  /* tensorflow::Tensor is always row-major */
              Xt.flat<Dtype>().data(),  /* ptr to data */
              X_dims[0],           /* num_rows */
              X_dims[1]           /* num_cols */);  //input_tensor.dim_size(1)

      auto mapOut = Eigen::Map<Eigen::Matrix<Dtype, //std::complex<double>
              Eigen::Dynamic,
              Eigen::Dynamic,
                Eigen::RowMajor>>(
              Zt->flat<Dtype>().data(),
                      X_dims[0],
                      X_dims[1] );

      //compute output  = exp[M*(-idt)] where M is the input matrix
      mapOut = (mapIn*Dtype(0,-delta_t_)).exp();

  }
};

//template struct MatrixExpFunctor<CPUDevice, int32>;
//template struct MatrixExpFunctor<CPUDevice, uint32>;
//template struct MatrixExpFunctor<CPUDevice, complex64>;
template struct MatrixExpFunctor<CPUDevice, std::complex<double>>;
template struct MatrixExpFunctor<CPUDevice, std::complex<float>>;


template <typename Dtype>
struct MatrixExpGrad<CPUDevice, Dtype> {
  static void launch(::tensorflow::OpKernelContext* ctx, const Tensor& topdiff_, const Tensor& Xt,
                     Tensor* grad_mA_, float delta_t_) {

      auto X_dims = (Xt.shape()).dim_sizes();

      const auto mapIn = Eigen::Map<const Eigen::Matrix<
      Dtype,           /* scalar element type */
      Eigen::Dynamic,  /* num_rows is a run-time value */
      Eigen::Dynamic,  /* num_cols is a run-time value */
      Eigen::RowMajor  >>(  /* tensorflow::Tensor is always row-major */
              Xt.flat<Dtype>().data(),  /* ptr to data */
              X_dims[0],           /* num_rows */
              X_dims[1]           /* num_cols */);  //input_tensor.dim_size(1)

      auto mapOut = Eigen::Map<Eigen::Matrix<Dtype, //std::complex<double>
              Eigen::Dynamic,
              Eigen::Dynamic,
      Eigen::RowMajor>>(
              grad_mA_->flat<Dtype>().data(),
                      X_dims[0],
                      X_dims[1] );

      const auto mapGradA = Eigen::Map<const Eigen::Matrix<
      Dtype,           /* scalar element type */
      Eigen::Dynamic,  /* num_rows is a run-time value */
      Eigen::Dynamic,  /* num_cols is a run-time value */
      Eigen::RowMajor  >>(  /* tensorflow::Tensor is always row-major */
              topdiff_.flat<Dtype>().data(),  /* ptr to data */
              X_dims[0],           /* num_rows */
              X_dims[1]           /* num_cols */);




      // the gradient with respect to the input matrix M is grad*A.adjoint()*(-idt) where A = exp[M*(-idt)]
      mapOut = mapGradA*(((mapIn*Dtype(0,-delta_t_)).exp()).adjoint())*Dtype(0,-delta_t_);

  }
};

//template struct MatrixExpGrad<CPUDevice, complex64>;
template struct MatrixExpGrad<CPUDevice, std::complex<double>>;
template struct MatrixExpGrad<CPUDevice, std::complex<float>>;


}  // namespace functor
}  // namespace tensorflow