// 2018, Cecilia Curreli

#include "tensorflow/core/framework/op.h"
#include "tensorflow/core/framework/op_kernel.h"

#include "matrix_exp_op.h"

#include <complex>


namespace tensorflow {

// Forward-Pass (CPU, GPU)
// --------------------------------------------------
template <typename Device, typename Dtype>
class MatrixExpOp : public OpKernel {
 public:
  explicit MatrixExpOp(OpKernelConstruction* ctx) : OpKernel(ctx) {
     OP_REQUIRES_OK(ctx, ctx->GetAttr("delta_t", &delta_t_));
     OP_REQUIRES(ctx, delta_t_ >= 0, tensorflow::errors::InvalidArgument("Need delta_t_ >= 0, got ", delta_t_ ));
  }

  void Compute(OpKernelContext* ctx) override {

      const tensorflow::Tensor& A = ctx->input(0);

      if (!ctx->status().ok()) {return;}

      TensorShape A_shape = A.shape();

      //check dimensions of input matrix
      OP_REQUIRES(ctx, tensorflow::TensorShapeUtils::IsSquareMatrix(A_shape),
                  tensorflow::errors::InvalidArgument("Input expects a 2-D matrix, but Input has Shape", A_shape));

      // Create an output tensor
      tensorflow::Tensor* Out = nullptr;
      OP_REQUIRES_OK(ctx, ctx->allocate_output(0, A_shape, &Out));

    ::tensorflow::functor::MatrixExpFunctor<Device, Dtype>::launch(ctx, A,  Out, delta_t_);
  }

 private:
  TF_DISALLOW_COPY_AND_ASSIGN(MatrixExpOp);
  float delta_t_;
};

// Backward-Pass (CPU, GPU)
// --------------------------------------------------
template <typename Device, typename Dtype>
class MatrixExpGradOp : public OpKernel {
 public:
  explicit MatrixExpGradOp(OpKernelConstruction* ctx) : OpKernel(ctx) {
      OP_REQUIRES_OK(ctx, ctx->GetAttr("delta_t", &delta_t_));
      OP_REQUIRES(ctx, delta_t_ >= 0, tensorflow::errors::InvalidArgument("Need delta_t_ >= 0, got ", delta_t_ ));
  }

  void Compute(OpKernelContext* ctx) override {
    const Tensor& X = ctx->input(1);
    const Tensor& topdiff = ctx->input(0);

    if (!ctx->status().ok()) { return; }

    Tensor* grad_X = nullptr;
    OP_REQUIRES_OK(ctx, ctx->allocate_output(0, X.shape(), &grad_X));

    ::tensorflow::functor::MatrixExpGrad<Device, Dtype>::launch(ctx, topdiff, X, grad_X, delta_t_);
  }

private:
float delta_t_;
};


#define REGISTER_CUSTOM_OP(NAME, DEVICE, T)                       \
  REGISTER_KERNEL_BUILDER(                                        \
      Name(#NAME).Device(DEVICE_##DEVICE).TypeConstraint<T>("T"), \
      NAME##Op<DEVICE##Device, T>)

REGISTER_CUSTOM_OP(MatrixExp, CPU, std::complex<double>);
REGISTER_CUSTOM_OP(MatrixExp, CPU, std::complex<float>);

REGISTER_CUSTOM_OP(MatrixExpGrad, CPU, std::complex<double>);
REGISTER_CUSTOM_OP(MatrixExpGrad, CPU, std::complex<float>);



#ifdef GOOGLE_CUDA
//REGISTER_CUSTOM_OP(MatrixExp, GPU, complex64);
REGISTER_CUSTOM_OP(MatrixExp, GPU, complex128);
//REGISTER_CUSTOM_OP(MatrixExpGrad, GPU, complex64);
REGISTER_CUSTOM_OP(MatrixExpGrad, GPU, complex128);
#endif  // GOOGLE_CUDA
#undef REGISTER_CUSTOM_OP
}  // namespace tensorflow
