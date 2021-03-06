/* Copyright 2019 The Spin-Scenario Authors. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#pragma once
#include <kernel/spinsys/spin_system.h>
using namespace ssl::spinsys;
#include <nlopt.hpp>


#ifdef TENSORFLOW_ENABLED
#include "tensorflow/core/public/session.h"
#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/core/framework/tensor.h"
#include "tensorflow/core/util/events_writer.h"
#include "tensorflow/core/framework/graph.pb.h"
#include "tensorflow/core/graph/default_device.h"
#include "tensorflow/core/graph/graph_def_builder.h"
using namespace tensorflow;
using namespace tensorflow::ops;
#endif //TENSORFLOW_ENABLED

namespace ssl {
namespace oc {
struct spin_sys_superop {
  sp_cx_mat L0; // free ham.
  std::vector<sp_cx_mat> L0s; // broadband free hams.
  vec nominal_offset; // corresponding broadband offsets.
  vec profile; // corresponding broadband coefficients.
  std::vector<std::vector<double>> grad_bb;
  vec phi_dist;
  sp_cx_mat R; // relaxation op.
  rf_ham rf_ctrl; // rf hams.
};

struct state_traj {
  sp_cx_vec *forward;
  sp_cx_vec *backward;
  state_traj() {
    forward = nullptr;
    backward = nullptr;
  }
  state_traj(size_t nsteps) {
    forward = new sp_cx_vec[nsteps + 1];
    backward = new sp_cx_vec[nsteps + 1];
  }
  void free() {
	  delete []forward;
	  delete []backward;
  }
};

struct op_traj {
  sp_cx_mat *forward;
  sp_cx_mat *backward;
  op_traj() {
    forward = nullptr;
    backward = nullptr;
  }
  op_traj(size_t nsteps) {
    forward = new sp_cx_mat[nsteps + 1];
    backward = new sp_cx_mat[nsteps + 1];
  }
};
} /* namespace oc */
} /* namespace ssl */
