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

#include "coop_grape.h"
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;

namespace ssl {
namespace oc {

struct coop_var {
  sp_cx_vec rho_bar;
  sp_cx_vec rho_lambda;
  vector<sp_cx_vec> obsrv_state;
  //vector<string> obsrv_state_expr;
  vec beta;
  vec gamma;
  void alloc(spin_system *sys, sp_cx_vec targ_state) {
    std::map<string, sp_cx_vec> obsrv_state_map;
    obsrv_state_map = sys->cartesian_basis_states();
    std::map<string, sp_cx_vec>::iterator iter;
    for (iter = obsrv_state_map.begin(); iter != obsrv_state_map.end(); iter++) {
      obsrv_state.push_back(norm_state(iter->second));
      //obsrv_state_expr.push_back(iter->first);
    }
    size_t n = obsrv_state.size();
    beta = vec::Zero(n);
    gamma = vec::Zero(n);
    rho_bar = obsrv_state[0];
    rho_lambda = obsrv_state[0];
    rho_bar.setZero();
    rho_lambda.setZero();

    for (size_t i = 0; i < obsrv_state.size(); i++)
      gamma[i] = transfer_fidelity(targ_state, obsrv_state[i]);
  }
  void calc_beta() {
    for (size_t i = 0; i < obsrv_state.size(); i++)
      beta[i] = transfer_fidelity(rho_bar, obsrv_state[i]);
  }

  void calc_lambdaN(size_t ncoop) {
    rho_lambda.setZero();
    for (size_t i = 0; i < obsrv_state.size(); i++)
      rho_lambda += (gamma[i] - beta[i]) * obsrv_state[i];
    rho_lambda *= 2.0 / double(ncoop);
    rho_lambda.coeffRef(0) = 1; // NOTE.
  }

  double calc_obj() {
    vec val = gamma - beta;
    return 1 - val.cwiseAbs2().sum();
  }
};

coop_grape::coop_grape(spin_system &sys) : grape(sys) {
  is_coop_= true;
}

coop_grape::~coop_grape() {
}

void coop_grape::assign_pulse(const sol::table &t) {
  grape::assign_pulse(t);
  size_t ncoop = retrieve_table_size_t("ncoop", t);

  string pattern = "rand_spline";
  if (is_retrievable("init_pattern", t))
    pattern = retrieve_table_str("init_pattern", t);

  std::map<string, rf_pattern>::const_iterator iter;
  iter = g_rf_pattern.find(pattern);

  for (size_t i = 0; i < ncoop; i++) {
    shaped_rf *ptr = rf_->Clone();
    ptr->copy_config_table();
    ptr->set_name("ms-coop " + to_string(i + 1)+ "-" + to_string(ncoop));
    ptr->set_shape(rf_->get_channels_str(), iter->second, 1);
    coop_rf_.push_back(ptr);
    //coop_rf_[i]->switch_rf_mode("amp/phase");
    //coop_rf_[i]->plot();
  }
}

void coop_grape::assign_aux_var() {
  omp_set_num_threads(omp_core_num);
  coop_var_ = new coop_var;
  coop_var_->alloc(sys_, targ_state_);
  for (size_t i = 0; i < coop_rf_.size(); i++)
    coop_traj_.push_back(state_traj(rf_->get_steps()));
}

sol::object coop_grape::optimize(const sol::table &t, sol::this_state s) {
  assign_state(t); // shared with grape.
  assign_nlopt(t); // shared with grape.
  assign_pulse(t);
  assign_aux_var();

  vector<double> x;
  for (size_t i = 0; i < coop_rf_.size(); i++) {
    coop_rf_[i]->convert2(_ux_uy);
    vector<double> raw = coop_rf_[i]->clone_raw_data();
    x.insert(x.end(), raw.begin(), raw.end()); // column major.
    raw.clear();
  }

  // optimization processing.
  nlopt::opt opt(opt_.algo, x.size());
  opt.set_xtol_rel(opt_.tol_f);
  if (opt_.max_eval > 0)
    opt.set_maxeval(opt_.max_eval);
  if (opt_.max_time > 0)
    opt.set_maxtime(opt_.max_time);
  if (opt_.stopval > 0)
    opt.set_stopval(opt_.stopval);

  opt.set_max_objective(co_objfunc, this);

  // boundary limitation.
  //opt_amplitude_constraint(opt);

  double max_f = 0;
  nlopt::result result = opt.optimize(x, max_f);
  if (result == nlopt::SUCCESS)
    ssl_color_text("info", "pulse optimization succeed.\n");
  if (result == nlopt::FAILURE)
    ssl_color_text("err", "pulse optimization failed.\n");
  if (result == nlopt::MAXEVAL_REACHED)
    ssl_color_text("warn", "pulse optimization terminated due to maximum iterations limitation.\n");

  opt_.max_f = max_f;

  sol::state_view lua(s);
  sol::table rfs = lua.create_table();
  int dim = (int) rf_->get_dims();
  for (size_t i = 0; i < coop_rf_.size(); i++) {
    coop_rf_[i]->update_raw_data(x.data() + i * dim);
    rfs.add((seq_block *) coop_rf_[i]);
  }
  h5write();
  return rfs;
}

void coop_grape::h5write(string file_name) const {
  if (file_name.empty()) {
    string time_s = sys_time();
    file_name = "ms_coop_" + time_s + ".h5";
  }
  H5File file(file_name, H5F_ACC_TRUNC);

  for (size_t i = 0; i < coop_rf_.size(); i++) {
    coop_rf_[i]->switch_rf_mode("amp/phase");
    coop_rf_[i]->h5write(file, coop_rf_[i]->name());
  }

  ssl::utility::h5write(file, nullptr, "obj", stl2vec(opt_.vf));
  file.close();
}

double coop_grape::co_objfunc(const vector<double> &x, vector<double> &g, void *func) {
  return ((coop_grape *) func)->co_objfunc(x, g);
}

double coop_grape::co_objfunc(const vector<double> &x, vector<double> &g) {
  int ncoop = (int) coop_rf_.size();
  int dim = (int) rf_->get_dims();
  Eigen::Map<const mat> shape(x.data(), dim, ncoop);
  Eigen::Map<mat> grad(g.data(), dim, ncoop);

  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  double dt = rf_->width_in_ms() * 1e-3 / double(nsteps); // into s.

  sp_cx_mat L0 = superop_.L0 + ci * superop_.R;

#pragma omp parallel for
  for (int p = 0; p < ncoop; p++) {
    coop_traj_[p].forward[0] = init_state_;

    sp_cx_mat L;
    double kx = 1, ky = 1;
    sp_cx_vec rho = coop_traj_[p].forward[0];
    for (size_t i = 0; i < nsteps; i++) {
      L = L0;
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(shape, p, i, j, nchannels, kx, ky);
      coop_traj_[p].forward[i + 1] = ssl::spinsys::step(rho, L, dt);
      rho = coop_traj_[p].forward[i + 1];
    }
  }

  coop_var_->rho_bar.setZero();
  for (int p = 0; p < ncoop; p++)
    coop_var_->rho_bar += coop_traj_[p].forward[nsteps];
  coop_var_->rho_bar /= (double) ncoop;

  coop_var_->rho_bar.coeffRef(0) = 0;
  coop_var_->rho_bar = norm_state(coop_var_->rho_bar);

  coop_var_->calc_beta();
  coop_var_->calc_lambdaN(ncoop);

#pragma omp parallel for
  for (int p = 0; p < ncoop; p++) {
    coop_traj_[p].backward[nsteps] = coop_var_->rho_lambda;

    sp_cx_mat L;
    double kx = 1, ky = 1;
    sp_cx_vec rho = coop_traj_[p].backward[nsteps];
    for (int i = nsteps - 1; i >= 0; i--) {
      L = L0.adjoint();
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(shape, p, i, j, nchannels, kx, ky);
      coop_traj_[p].backward[i] = ssl::spinsys::step(rho, L, -dt);
      rho = coop_traj_[p].backward[i];
    }
  }

#pragma omp parallel for
  for (int p = 0; p < ncoop; p++) {
    sp_cx_mat Gx, Gy, tmp;
    sp_cx_mat L;
    int k = 0;
    double kx = 1, ky = 1;
    for (size_t i = 0; i < nsteps; i++) {
      L = L0;
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(shape, p, i, j, nchannels, kx, ky);
      tmp = coop_traj_[p].backward[i + 1].adjoint() * ssl::spinsys::propagator(L, dt);
      for (size_t j = 0; j < nchannels; j++) {
        Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
        Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);
        grad(k, p) = traced(tmp * Gx * coop_traj_[p].forward[i]).real();
        grad(k + 1, p) = traced(tmp * Gy * coop_traj_[p].forward[i]).real();
        k += 2;
      }
    }
  }

  double val = coop_var_->calc_obj();
  std::cout << boost::format("==> %04d [%.7f]\n") % (++opt_.iter_count) % val;
  opt_.vf.push_back(val);
  return val;
}

void coop_grape::projection(const sol::table &t) {
  sol::object val = retrieve_table("init_state", t);
  sp_cx_vec init_state = sys_->smart_state(val.as<string>());
  init_state = norm_state(init_state);

  val = retrieve_table("coop_rf", t);
  sol::table rf_table = val.as<sol::table>();

  vector<shaped_rf *> co_rfs;
  for (size_t i = 0; i < rf_table.size(); i++) {
    sol::object val = rf_table[i + 1];
    seq_block *p_rf = val.as<seq_block *>();
    co_rfs.push_back((shaped_rf *) p_rf);
    //co_rfs[i]->switch_rf_mode("ux/uy");
  }

  size_t nsteps = co_rfs[0]->get_steps();
  size_t nchannels = co_rfs[0]->get_channels();
  double dt = co_rfs[0]->width_in_ms() * 1e-3 / double(nsteps); // into s.

  std::map<string, sp_cx_vec> obsrv_state_map;
  vector<string> expr;
  vector<sp_cx_vec> obsrv_state;
  if (is_retrievable("observ_states", t)) {
    val = retrieve_table("observ_states", t);
    sol::table expr_table = val.as<sol::table>();
    if (expr_table.size() == 0) {
      obsrv_state_map = sys_->cartesian_basis_states();
    } else {
      for (size_t i = 0; i < expr_table.size(); i++) {
        sol::object val = expr_table[i + 1];
        string exp = val.as<string>();
        sp_cx_vec rho = sys_->smart_state(exp);
        obsrv_state_map.insert(pair<string, sp_cx_vec>(exp, rho));
      }
    }
  } else {
    obsrv_state_map = sys_->cartesian_basis_states();
  }

  std::map<string, sp_cx_vec>::iterator iter;
  for (iter = obsrv_state_map.begin(); iter != obsrv_state_map.end(); iter++) {
    obsrv_state.push_back(norm_state(iter->second));
    expr.push_back(iter->first);
  }

  //for (size_t p = 0; p < expr.size(); p++)
    //std::cout << "row " << p + 1 << "=>" << expr[p] << "\n";

  string str_opt = "step";
  if (is_retrievable("option", t))
    str_opt = retrieve_table_str("option", t);

  int nrows = 0;
  if (str_opt == "step")
    nrows = nsteps;
  else if (str_opt == "broadband")
    nrows = superop_.L0s.size();
  else {
    string s = "unknown projection option ** " + str_opt + " ** using ‘step’ or 'broadband' instead.";
    throw std::runtime_error(s.c_str());
  }

  mat transfer_wave = mat::Zero(nrows, expr.size());

  if (str_opt == "step") {

    vector<state_traj> co_traj_rho;
    for (size_t i = 0; i < co_rfs.size(); i++)
      co_traj_rho.push_back(state_traj(nsteps));

    sp_cx_mat L0 = superop_.L0 + superop_.R;
    int ncoop = co_rfs.size();

    vector<double> x;
    for (size_t i = 0; i < co_rfs.size(); i++) {
      co_rfs[i]->convert2(_ux_uy);
      vector<double> raw = co_rfs[i]->clone_raw_data();
      x.insert(x.end(), raw.begin(), raw.end()); // column major.
      raw.clear();
    }

    Eigen::Map<const mat> shape(x.data(), co_rfs[0]->get_dims(), ncoop);

    omp_set_num_threads(omp_core_num);
#pragma omp parallel for
    for (int p = 0; p < ncoop; p++) {
      co_traj_rho[p].forward[0] = init_state;
      sp_cx_mat L;
      double kx = 1, ky = 1;
      sp_cx_vec rho = co_traj_rho[p].forward[0];
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(shape, p, i, j, nchannels, kx, ky);
        co_traj_rho[p].forward[i + 1] = ssl::spinsys::step(rho, L, dt);
        rho = co_traj_rho[p].forward[i + 1];
      }
    }

    state_traj rho_bar(nsteps);
    rho_bar.forward[0] = init_state;
#pragma omp parallel for
    for (int i = 0; i < (int) nsteps; i++) {
      rho_bar.forward[i + 1] = init_state;
      rho_bar.forward[i + 1].setZero();
      for (int p = 0; p < ncoop; p++)
        rho_bar.forward[i + 1] += co_traj_rho[p].forward[i + 1];

      rho_bar.forward[i + 1].coeffRef(0) = 0;
      rho_bar.forward[i + 1] = norm_state(rho_bar.forward[i + 1]);

      for (size_t k = 0; k < obsrv_state.size(); k++) {
        sp_cx_vec compo = obsrv_state[k];
        transfer_wave(i, k) = transfer_fidelity(rho_bar.forward[i + 1], compo) / transfer_fidelity(compo, compo);
      }
    }
  } else if (str_opt == "broadband") {
//        omp_set_num_threads(omp_core_num);

//        vector<sp_cx_vec*> forward_trajs_parfor = vector<sp_cx_vec*>(omp_core_num);
//        for (int i = 0; i < omp_core_num; i++) {
//            sp_cx_vec* forward = new sp_cx_vec[nsteps + 1];
//            forward[0] = init_state_;
//            forward_trajs_parfor[i] = forward;
//        }

//        #pragma omp parallel for
//        for (int p = 0; p < nrows; p++) {
//            int id = omp_get_thread_num();

//            sp_cx_mat L;
//            sp_cx_mat L0 = superop_.L0s[p] + ci*superop_.R;
//            sp_cx_vec rho = forward_trajs_parfor[id][0];

//            for (size_t i = 0; i < nsteps; i++) {
//                L = L0;
//                for (size_t j = 0; j < nchannels; j++)
//                    L += update_rf_ham(x.data(), i, j, nchannels);
//                forward_trajs_parfor[id][i + 1] = ssl::spinsys::step(rho, L, dt);
//                rho = forward_trajs_parfor[id][i + 1];
//            }

//            std::map<string, sp_cx_vec>::iterator iter;
//            for (size_t k = 0; k < expr.size(); k++) {
//                iter = obsrv_state_map.find(expr[k]);
//                if (iter != obsrv_state_map.end()) {
//                    sp_cx_vec compo = iter->second;
//                    transfer_wave(p, k) = transfer_fidelity(rho, compo) / transfer_fidelity(compo, compo);
//                }
//            }
//        }
  }

  // plot.
  sol::table lines = g_lua->create_table();
  for (int i = 0; i < transfer_wave.cols(); i++) {
    vec line = transfer_wave.col(i);
    lines.add(line);
  }

  //vec ss1 = transfer_wave.col(0).cwiseAbs2();
  //vec ss2 = transfer_wave.col(1).cwiseAbs2();
  //vec ss = ss1 + ss2;
  //vec s = ss.cwiseSqrt();
  //lines.add(s);

  string fig_spec;
  vec time;
  if (str_opt == "step") {
    time = vec::LinSpaced(nrows, 0, co_rfs[0]->width_in_ms());
    fig_spec = "title<transfer coefficient> xlabel<pulse duration / ms>";
  } else if (str_opt == "broadband") {
    int n = superop_.nominal_offset.size();
    time = vec::LinSpaced(nrows, superop_.nominal_offset[0], superop_.nominal_offset[n - 1]);
    fig_spec = "title<transfer coefficient> xlabel<freq offset / kHz>";
  }

  if (expr.size() > 5)
    fig_spec += " gnuplot<set key outside>";

  fig_spec += " color<PairedLines,16>";
  fig_spec += " lw<4>";

  string lege;
  for (size_t i = 0; i < expr.size(); i++)
    lege += expr[i] + ";";
  fig_spec += " legend<" + lege + ">";

  plot(fig_spec, line_series(time, lines));

  ofstream ofstr("basis_projection.dat");
  ofstr.precision(3);
  ofstr << transfer_wave;
  ofstr.close();

  return;

  vec2 xrange, yrange;
  if (str_opt == "step") {
    xrange[0] = 0;
    xrange[1] = co_rfs[0]->width_in_ms();
  } else if (str_opt == "broadband") {
    int n = superop_.nominal_offset.size();
    xrange[0] = superop_.nominal_offset[0];
    xrange[1] = superop_.nominal_offset[n - 1];
  }
  yrange[0] = 1;
  yrange[1] = expr.size();
  utility::map transfer_map(transfer_wave.transpose().matrix());
  transfer_map.xrange = xrange;
  transfer_map.yrange = yrange;
  (*g_lua)["_map"] = transfer_map;
  g_lua->script("plot('title<transfer coefficient> xlabel<pulse duration / ms> ylabel<projection component>', _map)");
}

sp_cx_mat coop_grape::update_rf_ham(Eigen::Map<const mat> &m,
                                    int scan,
                                    size_t step,
                                    size_t channel,
                                    size_t nchannels,
                                    double kx,
                                    double ky) const {
  double ux = m(2 * nchannels * step + 2 * channel, scan);
  double uy = m(2 * nchannels * step + 2 * channel + 1, scan);
  // Rf inhom
  ux *= kx;
  uy *= ky;
  return ux * superop_.rf_ctrl.Lx[channel] + uy * superop_.rf_ctrl.Ly[channel];
}
}
}
