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

#ifndef SPIN_SCENARIO_TF_OPT_H
#define SPIN_SCENARIO_TF_OPT_H

#include "oc_auxiliary.h"
#include <kernel/seq/shaped_rf.h>
using namespace ssl::seq;

namespace ssl {
namespace oc {

#ifdef TENSORFLOW_ENABLED
class tf_opt {

 protected:

 public:
  tf_opt(spin_system &sys);
  virtual ~tf_opt();
  seq_block &optimize(const sol::table &t);
  virtual void projection(const sol::table &t);
  double maxf() const;

 protected:

  void create_graph(std::string graph_type);
  bool call_system(int size, int n_channels, double dt, std::string graph_definition, std::string graph_type);
  void initialize_placeholder();
  void create_tf_variables_list();
  void update_tf_variables(const double *x);

  static double objfunc(const std::vector<double> &x, std::vector<double> &grad, void *func);
  double objfunc(const std::vector<double> &x, std::vector<double> &grad);

  void assign_state(const sol::table &t);
  void assign_nlopt(const sol::table &t);
  virtual void assign_pulse(const sol::table &t);
  virtual void assign_aux_var();
  virtual void h5write(std::string file_name = "") const;
//            void opt_amplitude_constraint(nlopt::opt &opt);
//            static double amplitude_constraint(unsigned n, const double *x, double *grad, void *data);
  void cartesian2polar(double amp, double phase, double gx, double gy, double &g_amp, double &g_phase);
  sp_cx_mat update_rf_ham(const double *x,
                          size_t step,
                          size_t channel,
                          size_t nchannels,
                          double kx = 1,
                          double ky = 1) const;
 protected:

  size_t nsteps;
  double dt; //in s
  bool relax;

  spin_sys_superop superop_;
  shaped_rf *rf_;

  sp_cx_vec init_state_;
  sp_cx_vec targ_state_;
  std::vector<sp_cx_vec> targ_list_;

  state_traj traj_;
  std::vector<state_traj> traj_omp_;

  opt_ctrl opt_;
  spin_system *sys_;

  Session *session; //std::unique_ptr<Session> session(NewSession(options)); //Yan
  GraphDef graph_def;
  SessionOptions opts;
  std::vector<Tensor> outputs;

  std::vector<std::pair<std::string, Tensor> > variables_list;
  std::vector<std::pair<std::string, Tensor> > placeholder_list;
  std::vector<std::string> gradient_names;
};
#endif  // TENSORFLOW_ENABLED

}
}
#endif //SPIN_SCENARIO_TF_OPT_H
