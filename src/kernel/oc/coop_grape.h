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

struct coop_var;
class coop_grape : public grape {
 public:
  coop_grape(spin_system &sys);
  virtual ~coop_grape();
  sol::object optimize(const sol::table &t, sol::this_state s);
  virtual void projection(const sol::table &t);
 protected:
  virtual void h5write(string file_name = "") const;
  virtual void assign_pulse(const sol::table &t);
  virtual void assign_aux_var();
  static double co_objfunc(const vector<double> &x, vector<double> &g, void *func);
  double co_objfunc(const vector<double> &x, vector<double> &g);

  sp_cx_mat update_rf_ham(Eigen::Map<const mat> &m,
                          int scan,
                          size_t step,
                          size_t channel,
                          size_t nchannels,
                          double kx = 1,
                          double ky = 1) const;
 private:
  vector<shaped_rf *> coop_rf_;
  vector<state_traj> coop_traj_;
  coop_var *coop_var_;
};
}
}
#endif // COOP_GRAPE_H
