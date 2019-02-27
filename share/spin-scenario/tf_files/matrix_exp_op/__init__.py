# 2018, Cecilia Curreli

# manually generated file
import tensorflow as tf
import os
from tensorflow.python.framework import ops

__all__ = ['matrix_exp']

path = os.path.join(os.path.dirname(__file__), 'matrix_exp_op.so')
_matrix_exp_module = tf.load_op_library(path)

matrix_exp = _matrix_exp_module.matrix_exp
_matrix_exp_grad = _matrix_exp_module.matrix_exp_grad


@ops.RegisterGradient("MatrixExp")
def _MatrixAddGrad(op, *grads):
  delta_t = op.get_attr('delta_t')
  matA = op.inputs[0]
  #matB = op.inputs[1]
  # top = op.outputs[0]
  topdiff = grads[0]
  return _matrix_exp_grad(topdiff, matA, delta_t=delta_t)
