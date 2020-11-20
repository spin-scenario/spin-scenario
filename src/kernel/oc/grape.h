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
extern int g_oc_iteration_save;
void set_oc_iteration_save(const sol::table &t);
enum opt_model {_rho2rho, _propagator};
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
  virtual void assign_x();
  virtual void assign_aux_var();
  virtual void print() const;
  virtual void h5write(std::string file_name = "") const;
  static double objfunc_broadband(const std::vector<double> &x, std::vector<double> &grad, void *func);
  static double objfunc_propagator(const std::vector<double> &x, std::vector<double> &grad, void *func);
  // double objfunc(const std::vector<double> &x, std::vector<double> &grad);
  double objfunc_broadband(const std::vector<double> &x, std::vector<double> &grad);
  double objfunc_propagator(const std::vector<double> &x, std::vector<double> &grad);
  void iteration_shape(const std::vector<double> &x, int iter);

   sp_cx_mat update_rf_ham(const double *x,
                          size_t step,
                          size_t channel,
                          std::string ch_str,
                          size_t nchannels,
                          double kx = 1,
                          double ky = 1);
  sp_cx_mat propagator_derivative(const sp_cx_mat &L, const sp_cx_mat &Lrf, double dt);

  void assign_constraint(const sol::table &t);
  static double amplitude_constraint(unsigned n, const double *x, double *grad, void *data);

  //void cartesian2polar(double amp, double phase, double gx, double gy, double &g_amp, double &g_phase);

 protected:
  nlopt::opt *optimizer_;
  std::vector<double> x_;
  size_t x_dim_;

  limit_axis axis_;
  double max_val_;
  std::vector<double> obj_val_;



  opt_model opt_model_;
  spin_sys_superop superop_;
  shaped_rf *rf_;

  sp_cx_vec init_state_;
  sp_cx_vec targ_state_;
  std::vector<sp_cx_vec> targ_list_;

  //state_traj traj_;
  std::vector<state_traj> traj_omp_;

  sp_cx_mat init_op_;
  sp_cx_mat targ_op_;
  std::vector<sp_cx_mat> targ_op_list_;

  std::vector<op_traj> op_traj_omp_;

  spin_system *sys_;
  double epsilon_;   // rf scaling factor. [-epsilon_, epsilon_]
  int epsilon_num_;  // rf scaling number.
};

} /* namespace oc */
} /* namespace ssl */
