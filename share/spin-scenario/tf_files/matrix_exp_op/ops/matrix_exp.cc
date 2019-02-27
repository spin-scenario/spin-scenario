// 2018, Cecilia Curreli

#include "tensorflow/core/framework/op.h"
#include "tensorflow/core/framework/shape_inference.h"

namespace tensorflow {

using ::tensorflow::shape_inference::InferenceContext;
using ::tensorflow::shape_inference::ShapeHandle;

namespace shape_inference {
Status UnchangedShape(InferenceContext* c) {
  c->set_output(0, c->input(0));
  return Status::OK();
}
}  // namespace shape_inference

REGISTER_OP("MatrixExp")
    .Attr("delta_t: float")
    .Input("matrix_a: T")
    .Output("output_matrix: T")
    .Attr("T: {complex128, complex64}")
    .SetShapeFn([](::tensorflow::shape_inference::InferenceContext* c) {

    ::tensorflow::shape_inference::ShapeHandle input_shape;
    TF_RETURN_IF_ERROR(c->WithRank(c->input(0), 2, &input_shape)); //If <shape> has rank <rank> return OK
    c->set_output(0, c->input(0));

    return tensorflow::Status::OK();
})
    .Doc(R"doc(
Computes the matrix exponential of the input matrix multiplied by -i*delta_t.

matrix_a: matrix [M, M].
delta_t: An additional constant term.
output_matrix: a matrix of the same size [M, M] containing exp(-i*delta_t*A).

)doc");

REGISTER_OP("MatrixExpGrad")
    .Attr("delta_t: float")
    .Input("grad: T")
    .Input("matrix_a: T")
    .Output("grad_a: T")
    .Attr("T: {complex128, complex64} ")
    .SetShapeFn([](InferenceContext* c) {
      c->set_output(0, c->input(1)); // grad_aA has same shape as matrix_A
      return ::tensorflow::Status::OK();
    })
    .Doc(R"doc(
Returns gradient of exp(-i*delta_t*A) with respect to A.
)doc");

//.Attr("T: complexnumbertype = DT_COMPLEX128")

}  // namespace tensorflow
