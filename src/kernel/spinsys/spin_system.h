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

#include "spin_system_auxiliary.h"

namespace ssl {
namespace spinsys {

class spin_system {
 public:
  spin_system(const std::string &filename = "");
  spin_system(const sol::table &t);
  ~spin_system();

  void set_sys(const sol::table &t);
  void set_isotopes(const std::string list); // "1H, 13C"
  void set_magnet_field(double tesla);
  void set_proton_freq(double MHz);
  double get_proton_freq() const;

  inline double get_magnet_field() const {
    return comp_.get_field();
  }
  std::string isotopes() const {
    return comp_.isotopes();
  }
  std::vector<std::string> channels() const {
    return comp_.channels();
  }

  inline int nspins() const {
    return comp_.nspins();
  }

  sp_mat op(const std::string &list, op_side type) const;
  sp_mat op(const sol::table &t, op_side type) const;
  sp_mat op(const std::string &list) const;
  sp_mat op(const sol::table &t) const;
  sp_cx_mat smart_op(const std::string expr) const;

  sp_cx_vec state(const std::string &list) const;
  sp_cx_vec state(const sol::table &t) const;
  sp_cx_vec smart_state(const std::string expr) const;

  sp_cx_vec equilibrium_state() const;

  vec nominal_broadband();

  std::map<std::string, sp_cx_vec> cartesian_basis_states() const;
  ham_op hamiltonian(const op_side type, bool build_aniso) const;
  sp_cx_mat relaxation() const;
  sp_cx_mat free_hamiltonian() const;
  sp_cx_mat total_hamiltonian() const;
  rf_ham rf_hamiltonian() const;
  std::vector<sp_cx_mat> free_hamiltonians();
  sol::object free_hamiltonians(sol::this_state s);
 private:
  composition comp_;  // determines the composition of the spin system.
  interaction inter_;
};

}
}
