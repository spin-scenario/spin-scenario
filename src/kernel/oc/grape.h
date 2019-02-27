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
#include "oc_auxiliary.h"
#include <kernel/seq/shaped_rf.h>
using namespace ssl::seq;

namespace ssl {
namespace oc {
class grape {
 public:
  grape(spin_system &sys);
  virtual ~grape();
  seq_block &optimize(const sol::table &t);
  virtual void projection(const sol::table &t);
  double maxf() const;
 protected:
  void assign_state(const sol::table &t);
  void assign_nlopt(const sol::table &t);
  virtual void assign_pulse(const sol::table &t);
  virtual void assign_aux_var();
  virtual void h5write(string file_name = "") const;
  static double objfunc(const vector<double> &x, vector<double> &grad, void *func);
  static double objfunc_broadband(const vector<double> &x, vector<double> &grad, void *func);
  double objfunc(const vector<double> &x, vector<double> &grad);
  double objfunc_broadband(const vector<double> &x, vector<double> &grad);

  sp_cx_mat update_rf_ham(const double *x,
                          size_t step,
                          size_t channel,
                          size_t nchannels,
                          double kx = 1,
                          double ky = 1) const;
  sp_cx_mat propagator_derivative(const sp_cx_mat &L, const sp_cx_mat &Lrf, double dt);

  void opt_amplitude_constraint(nlopt::opt &opt);
  static double amplitude_constraint(unsigned n, const double *x, double *grad, void *data);

  void cartesian2polar(double amp, double phase, double gx, double gy, double &g_amp, double &g_phase);
 protected:
  spin_sys_superop superop_;
  shaped_rf *rf_;

  sp_cx_vec init_state_;
  sp_cx_vec targ_state_;
  vector<sp_cx_vec> targ_list_;

  state_traj traj_;
  vector<state_traj> traj_omp_;

  opt_ctrl opt_;
  spin_system *sys_;
  bool is_coop_;
};

} /* namespace oc */
} /* namespace ssl */
