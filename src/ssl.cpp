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

#include "ssl.h"

namespace ssl {

void bindings(sol::state &lua) {
  // ssl IO utility bindings.
  sol::table ssl = lua.create_named_table("ssl");

  lua.set_function("zero_padding", &ssl::utility::set_zero_pad_api);
  lua.set_function("apodization", &ssl::utility::set_decay_rate_api);
  lua.set_function("get", &ssl::utility::get);
  lua.set_function("max", &ssl::utility::max);
  lua.set_function("cas", &ssl::utility::yacas_evaluate);
  lua.set_function("core", &ssl::utility::set_openmp_core);
 
  ssl.set_function(
      "print",
      sol::overload(
          sol::resolve<void(sol::variadic_args, const sp_mat &)>(print<double>),
          sol::resolve<void(sol::variadic_args, const sp_cx_mat &)>(print<cd>),
          sol::resolve<void(sol::variadic_args, const sp_vec &)>(print<double>),
          sol::resolve<void(sol::variadic_args, const sp_cx_vec &)>(print<cd>),
          sol::resolve<void(sol::variadic_args, const mat &)>(print<double>),
          sol::resolve<void(sol::variadic_args, const cx_mat &)>(print<cd>),
          sol::resolve<void(sol::variadic_args, const vec &)>(print<double>),
          sol::resolve<void(sol::variadic_args, const ivec &)>(print<int>),
          sol::resolve<void(sol::variadic_args, const cx_vec &)>(print<cd>),
          sol::resolve<void(sol::variadic_args, const seq_block &)>(
              ssl::seq::print)));

  ssl.set_function(
      "write",
      sol::overload(
          sol::resolve<void(std::string, sol::variadic_args, const sp_mat &)>(write),
          sol::resolve<void(std::string, sol::variadic_args, const sp_cx_mat &)>(
              write),
          sol::resolve<void(std::string, sol::variadic_args, const sp_vec &)>(write),
          sol::resolve<void(std::string, sol::variadic_args, const sp_cx_vec &)>(
              write),
          sol::resolve<void(std::string, sol::variadic_args, const mat &)>(write),
          sol::resolve<void(std::string, sol::variadic_args, const cx_mat &)>(write),
          sol::resolve<void(std::string, sol::variadic_args, const vec &)>(write),
          sol::resolve<void(std::string, sol::variadic_args, const cx_vec &)>(write),
          sol::resolve<void(std::string, sol::variadic_args, const seq_block &)>(
              ssl::seq::write)));

  ssl.set_function(
      "specgram",
      sol::overload(
          sol::resolve<void(const sol::table &)>(&ssl::seq::specgram),
          sol::resolve<void(std::string, const sol::table &)>(&ssl::seq::specgram)));

  lua.set("ci", ci);

  lua.set_function("h5read_mat",  sol::overload(sol::resolve<mat(std::string, std::string)>(h5read_mat)));

  lua.set_function("table2mat", table2mat);
  lua.set_function("table2vec", table2vec);
  lua.set_function("vec2table", vec2table);

  lua.set_function("pw90", set_pw90_api);

  lua.set_function("peak_grad", set_max_grad_api);
  lua.set_function("slew_rate", set_max_slew_rate_api);
  // lua.set_function("spec_grad", set_grad);

  lua.set_function("B0", set_B0_api);

  lua.set_function("phase_cycle_steps", set_phase_cycle_steps_api);

  // lua["pw90"] = &pw90;

  lua.new_usertype<sp_mat>(
      "sp_mat", sol::meta_function::addition, &ssl::utility::add<sp_mat>,
      sol::meta_function::subtraction, &ssl::utility::sub<sp_mat>,
      sol::meta_function::multiplication,
      sol::overload(
          sol::resolve<sp_mat(const sp_mat &, const sp_mat &)>(
              &ssl::utility::mul<sp_mat>),
          sol::resolve<sp_mat(const sp_mat &, double)>(
              &ssl::utility::mul<sp_mat>),
          sol::resolve<sp_mat(double, const sp_mat &)>(
              &ssl::utility::mul<sp_mat>),
          sol::resolve<sp_cx_mat(const sp_mat &, cd)>(&ssl::utility::mul),
          sol::resolve<sp_cx_mat(cd, const sp_mat &)>(&ssl::utility::mul)),
      sol::meta_function::division,
      sol::overload(
          sol::resolve<sp_mat(const sp_mat &, double)>(
              &ssl::utility::div<sp_mat>),
          sol::resolve<sp_cx_mat(const sp_mat &, cd)>(&ssl::utility::div)));
  lua.new_usertype<sp_cx_mat>(
      "sp_cx_mat", sol::meta_function::addition, &ssl::utility::add<sp_cx_mat>,
      sol::meta_function::subtraction, &ssl::utility::sub<sp_cx_mat>,
      sol::meta_function::multiplication,
      sol::overload(
          sol::resolve<sp_cx_mat(const sp_cx_mat &, const sp_cx_mat &)>(
              &ssl::utility::mul<sp_cx_mat>),
          sol::resolve<sp_cx_mat(const sp_cx_mat &, double)>(
              &ssl::utility::mul<sp_cx_mat>),
          sol::resolve<sp_cx_mat(double, const sp_cx_mat &)>(
              &ssl::utility::mul<sp_cx_mat>),
          sol::resolve<sp_cx_mat(const sp_cx_mat &, cd)>(&ssl::utility::mul),
          sol::resolve<sp_cx_mat(cd, const sp_cx_mat &)>(&ssl::utility::mul)),
      sol::meta_function::division,
      sol::overload(
          sol::resolve<sp_cx_mat(const sp_cx_mat &, double)>(
              &ssl::utility::div<sp_cx_mat>),
          sol::resolve<sp_cx_mat(const sp_cx_mat &, cd)>(&ssl::utility::div)));
  lua.new_usertype<sp_cx_vec>(
      "sp_cx_vec", sol::meta_function::addition, &ssl::utility::add<sp_cx_vec>,
      sol::meta_function::subtraction, &ssl::utility::sub<sp_cx_vec>,
      sol::meta_function::multiplication,
      sol::overload(
          sol::resolve<sp_cx_vec(const sp_cx_vec &, double)>(
              &ssl::utility::mul<sp_cx_vec>),
          sol::resolve<sp_cx_vec(double, const sp_cx_vec &)>(
              &ssl::utility::mul<sp_cx_vec>),
          sol::resolve<sp_cx_vec(const sp_cx_vec &, cd)>(&ssl::utility::mul),
          sol::resolve<sp_cx_vec(cd, const sp_cx_vec &)>(&ssl::utility::mul)),
      sol::meta_function::division,
      sol::overload(
          sol::resolve<sp_cx_vec(const sp_cx_vec &, double)>(
              &ssl::utility::div<sp_cx_vec>),
          sol::resolve<sp_cx_vec(const sp_cx_vec &, cd)>(&ssl::utility::div)));

  lua.set_function(
      "fft_1d",
      sol::overload(
          sol::resolve<cx_vec(const cx_vec &)>(ssl::utility::fft_1d),
          sol::resolve<cx_vec(const sol::table &)>(ssl::utility::fft_1d)));

  // lua.set_function("linspace", &ssl::utility::linspace_lua);
  lua.set_function("linspace", &ssl::utility::linspace);
  lua.set_function("random", &ssl::utility::random);
  lua.set_function("mat", &ssl::utility::vec_table);

  ssl.new_usertype<line>(
      "line", sol::constructors<
                  line(std::string), line(std::string, std::string), line(const vec &),
                  line(const vec &, std::string), line(const vec &, const vec &),
                  line(const vec &, const vec &, std::string),
                  line(const sol::table &), line(const sol::table &, std::string),
                  line(const sol::table &, const sol::table &),
                  line(const sol::table &, const sol::table &, std::string),
                  line(const sol::table &, const sol::table &, std::string)>());
  ssl.new_usertype<line_series>(
      "line_series",
      sol::constructors<line_series(std::string), line_series(const sol::table &),
                        line_series(const vec &, const sol::table &)>());
  ssl.new_usertype<utility::map>(
      "map",
      sol::constructors<utility::map(std::string), utility::map(std::string, std::string),
                        utility::map(const mat &),
                        utility::map(const mat &, std::string)>());

  ssl.set_function(
      "plot",
      sol::overload(
          // sol::resolve<void(const sol::table&)>(ssl::seq::plot),
          sol::resolve<void(sol::variadic_args, const seq_block &)>(
              ssl::seq::plot),
          sol::resolve<void(sol::variadic_args, const line &)>(
              ssl::utility::plot),
          sol::resolve<void(std::string, sol::variadic_args, const line &)>(
              ssl::utility::plot),
          sol::resolve<void(std::string, const line_series &)>(
              ssl::utility::plot),
          sol::resolve<void(sol::variadic_args, const utility::map &)>(
              ssl::utility::plot),
          sol::resolve<void(std::string, sol::variadic_args, const utility::map &)>(
              ssl::utility::plot),
          sol::resolve<void(const line_series &)>(ssl::utility::plot)
          //sol::resolve<void(std::string, const line_series &)>(ssl::utility::plot)
      ));

    lua.set_function("output_terminal", set_output_terminal);


  // pulse sequence bindings.
  ssl.new_usertype<seq_block_factory>(
      "seq_block_factory", sol::constructors<sol::types<>>(), "create",
      sol::resolve<sol::object(const sol::table &, sol::this_state) const>(
          &seq_block_factory::clone_seq_block));

  // oc bindings.
#ifdef TENSORFLOW_ENABLED
  ssl.new_usertype<tf_opt>(
      "cecilia", sol::constructors<sol::types<const spin_system &>>(),
      "projection", &tf_opt::projection, "optimize", &tf_opt::optimize);
#endif  // TENSORFLOW_ENABLED
  
  ssl.new_usertype<grape>(
      "rf_optimizer", sol::constructors<sol::types<const spin_system &>>(),
      "maxf", sol::property(&grape::maxf), "projection", &grape::projection,
      "optimize",
      sol::overload(
          sol::resolve<seq_block &(const sol::table &)>(&grape::optimize)));

  ssl.new_usertype<coop_grape>(
      "multi_rf_optimizer",
      sol::constructors<sol::types<const spin_system &>>(),
      "maxf", sol::property(&grape::maxf),
      "projection", &coop_grape::projection, "optimize",
      sol::overload(
          sol::resolve<sol::object(const sol::table &, sol::this_state)>(
              &coop_grape::optimize)));

    lua.set_function("spec_avg", spec_avg);
    lua.set_function("abs", &ssl::oc::abs);

  g_lua = &lua;

  ssl.new_usertype<seq_block>(
      "seq_block", sol::constructors<sol::types<>>(), "config",
      &seq_block::config, "set_name", &seq_block::set_name, "name",
      sol::property(&seq_block::name), "tau",
      sol::property(&seq_block::width_in_ms), "switch",
      &seq_block::switch_rf_mode, "rf_power",
      sol::property(&seq_block::rf_power), sol::meta_function::addition,
      &ssl::seq::concurrent, sol::meta_function::multiplication,
      &ssl::seq::set_cycle_priority, sol::meta_function::division,
      &ssl::seq::set_loop_array, sol::meta_function::modulus,
      sol::overload(
          sol::resolve<seq_block &(seq_block &, double)>(&ssl::seq::set_align),
          sol::resolve<seq_block &(seq_block &, std::string)>(
              &ssl::seq::set_align)));

  lua.set_function("area", &ssl::seq::area);
  lua.set_function("seq",
                   sol::overload(sol::resolve<seq_block &(const sol::table &)>(
                       &ssl::seq::serial)));
  lua.set_function("run",
                   sol::overload(sol::resolve<sol::object(const sol::table &)>(
                       &ssl::seq::run_seq_api)));
  lua.set_function("computeEngine", &ssl::seq::init_compute_engine);
  lua.set_function("seq_parm", &ssl::seq::init_seq_param);
  lua.set_function("reduce_phantom", &ssl::utility::reduce_phantom);

  lua.set_function("multiRF", &ssl::seq::multi_shaped_rf);
  lua.set_function("set_ham_Jcoupl", &ssl::physx::set_ham_Jcoupl);
  

  // phantom bindings.
  ssl.new_usertype<engine>(
      "engine", sol::constructors<sol::types<const sol::table &>>()
      //, "dim", sol::property(&phantom::dim)
  );

  // phantom bindings.
  ssl.new_usertype<phantom>(
      "phantom", sol::constructors<sol::types<>, sol::types<const char *>>(),
      "view", sol::resolve<void(const sol::table &) const>(&phantom::view)
      //"dim", sol::property(&phantom::dim)
  );

  // spin_system bindings.
  ssl.new_usertype<spin_system>(
      "spin_system",
      sol::constructors<sol::types<const std::string &>,
                        sol::types<const sol::table &>>(),
      "isotopes", sol::property(&spin_system::isotopes), "hamiltonian",
      &spin_system::hamiltonian, "freq",
      sol::property(&spin_system::get_proton_freq),"equilibrium_state",
      sol::property(&spin_system::equilibrium_state), "L0",
      sol::property(&spin_system::free_hamiltonian), "R",
      sol::property(&spin_system::relaxation), "L",
      sol::property(&spin_system::total_hamiltonian), "L0s",
      sol::overload(sol::resolve<sol::object(sol::this_state)>(
          &spin_system::free_hamiltonians)),
      "expr_op", &spin_system::smart_op, "basis_states",
      &spin_system::cartesian_basis_states, "op",
      sol::overload(
          sol::resolve<sp_mat(const sol::table &) const>(&spin_system::op),
          sol::resolve<sp_mat(const sol::table &, op_side) const>(
              &spin_system::op),
          sol::resolve<sp_mat(const std::string &) const>(&spin_system::op),
          sol::resolve<sp_mat(const std::string &, op_side) const>(
              &spin_system::op)),
      "expr_state", &spin_system::smart_state, "state",
      sol::overload(
          sol::resolve<sp_cx_vec(const sol::table &) const>(
              &spin_system::state),
          sol::resolve<sp_cx_vec(const std::string &) const>(&spin_system::state)));

  // lua.set_function("seq", &ssl::seq::seq);

  // lua.set_function("step", &ssl::spinsys::step);
  lua.set_function("acquire", &ssl::spinsys::acquire);
  // lua.set_function("projection", &ssl::spinsys::projection);

  lua.set_function("commu", &ssl::spinsys::commutator);
  lua.set_function("normalized",
                   sol::overload(sol::resolve<sp_cx_vec(const sp_cx_vec &)>(
                       normalized<cd>)));
}
}

