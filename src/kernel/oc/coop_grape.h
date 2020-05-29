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

#ifndef COOP_GRAPE_H
#define COOP_GRAPE_H
#include "grape.h"

namespace ssl {
namespace oc {

vec spec_avg(const sol::table &t);
vec abs(const vec &a, const vec &b);

enum coop_model { _ms_coop, _ss_coop};

struct ms_coop_ctrl;
class coop_grape : public grape {
 public:
  coop_grape(spin_system &sys);
  virtual ~coop_grape();
  sol::object optimize(const sol::table &t, sol::this_state s);
  virtual void projection(const sol::table &t);
 protected:
  virtual void h5write(std::string file_name = "") const;
  void assign_constraint(const sol::table &t);
  virtual void assign_pulse(const sol::table &t);
  virtual void assign_x();
  virtual void assign_aux_var();
  static double objfunc_broadband(const std::vector<double> &x, std::vector<double> &g, void *func);
  //static double objfunc_propagator(const std::vector<double> &x, std::vector<double> &grad, void *func); // reserved for propagator optimization case.
  double objfunc_broadband_ms_coop(const std::vector<double> &x, std::vector<double> &g);
  double objfunc_broadband_ss_coop(const std::vector<double> &x, std::vector<double> &g);
  double co_objfunc(const std::vector<double> &x, std::vector<double> &g); // this is for superposition state preparation, should be merged with objfunc_broadband.

  void assign_insert_propagator(const sol::table &t);
  virtual void print() const;
  sp_cx_mat update_rf_ham(Eigen::Map<const mat> &m, int scan, size_t step, size_t channel, std::string ch_str, size_t nchannels, double kx = 1, double ky = 1);
  sp_cx_mat update_rf_ham(Eigen::Map<const vec> &v, size_t step, size_t channel, std::string ch_str, size_t nchannels, double kx = 1, double ky = 1);
 private:
  std::vector<shaped_rf*> coop_rf_;
  std::vector<state_traj> coop_traj_;


  std::vector<sp_cx_mat> inserted_ops_; // only for ss-coop.

  ms_coop_ctrl* coop_par_;
  static coop_model coop_model_;
};
}
}
#endif // COOP_GRAPE_H
