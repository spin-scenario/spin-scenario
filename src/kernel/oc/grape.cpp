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

#include "grape.h"
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

namespace ssl {
namespace oc {

typedef struct {
  size_t i; // step index.
  size_t j; // channel index.
  size_t chs; // channel number.
  double amp;
} amp_constraint_data;

grape::grape(spin_system &sys)
    : rf_(nullptr) {
  is_coop_ = false;
  sys_ = &sys;
  superop_.rf_ctrl = sys.rf_hamiltonian();
  superop_.L0 = sys.free_hamiltonian();
  superop_.R = sys.relaxation();
  superop_.L0s = sys.free_hamiltonians();
  int bb_size = superop_.L0s.size();
  if (bb_size) {
    superop_.profile = vec::Ones(bb_size);
    superop_.nominal_offset = sys.nominal_broadband();
    superop_.grad_bb = vector<vector<double>>(bb_size);
  }
}

grape::~grape() {
}

double grape::maxf() const {
  return opt_.max_f;
}

seq_block &grape::optimize(const sol::table &t) {
  assign_state(t);
  assign_nlopt(t);
  assign_pulse(t);
  assign_aux_var();

  rf_->convert2(_ux_uy);
  vector<double> x = rf_->clone_raw_data();

  // optimization processing.
  nlopt::opt opt(opt_.algo, rf_->get_dims());
  opt.set_xtol_rel(opt_.tol_f);
  if (opt_.max_eval > 0)
    opt.set_maxeval(opt_.max_eval);
  if (opt_.max_time > 0)
    opt.set_maxtime(opt_.max_time);
  if (opt_.stopval > 0)
    opt.set_stopval(opt_.stopval);

  if (superop_.L0s.size() > 1)
    opt.set_max_objective(objfunc_broadband, this);
  else
    opt.set_max_objective(objfunc, this);

  // boundary limitation.
  opt_amplitude_constraint(opt);

  double max_f;
  nlopt::result result = opt.optimize(x, max_f);
  if (result == nlopt::SUCCESS)
    ssl_color_text("info", "pulse optimization succeed.\n");
  if (result == nlopt::FAILURE)
    ssl_color_text("err", "pulse optimization failed.\n");
  if (result == nlopt::MAXEVAL_REACHED)
    ssl_color_text("warn", "pulse optimization terminated due to maximum iterations limitation.\n");

  opt_.max_f = max_f;
  rf_->update_raw_data(x.data());

  grape::h5write();
  return *rf_;
}

double grape::objfunc(const vector<double> &x, vector<double> &grad, void *func) {
  return ((grape *) func)->objfunc(x, grad);
}

double grape::objfunc_broadband(const vector<double> &x, vector<double> &grad, void *func) {
  return ((grape *) func)->objfunc_broadband(x, grad);
}

void grape::assign_state(const sol::table &t) {
  init_state_ = sys_->smart_state(retrieve_table_str("init_state", t));
  init_state_ = norm_state(init_state_);

  sol::object obj = retrieve_table("targ_state", t);
  if (obj.get_type() == sol::type::string) {
    targ_state_ = sys_->smart_state(retrieve_table_str("targ_state", t));
    targ_state_ = norm_state(targ_state_);

    // in case of unified targ for diff freq offsets.
    if (superop_.L0s.size())
      targ_list_ = vector<sp_cx_vec>(superop_.L0s.size(), targ_state_);
  }

 if (obj.get_type() == sol::type::table) {
   targ_list_.clear();
    sol::table targ_table = obj.as<sol::table>();
    for (size_t i = 0; i < targ_table.size(); i++) {
      sol::object val = targ_table[i + 1];
      sp_cx_vec targ = val.as<sp_cx_vec>();
      targ=norm_state(targ);
      targ_list_.push_back(targ);
    }
  }
}

void grape::assign_nlopt(const sol::table &t) {
  // default algorithm for single pulse.
  opt_.algo = nlopt::LD_LBFGS; // default algorithm.

  // default para for COOP.
  if (is_coop_) {
    opt_.algo = nlopt::LD_MMA;
    opt_.tol_f = 1e-7;
  }

  if (is_retrievable("algorithm", t)) {
    string oc = retrieve_table_str("algorithm", t);
    boost::to_upper(oc);

    if (oc == "LBFGS")
      opt_.algo = nlopt::LD_LBFGS;

    if (oc == "MNA")
      opt_.algo = nlopt::LD_MMA;
  }

  if (is_retrievable("max_eval", t))
    opt_.max_eval = retrieve_table_int("max_eval", t);

  if (is_retrievable("max_time", t))
    opt_.max_time = retrieve_table_double("max_time", t);

  if (is_retrievable("tol_f", t))
    opt_.tol_f = retrieve_table_double("tol_f", t);

  if (is_retrievable("stopval", t))
    opt_.stopval = retrieve_table_double("stopval", t);

  if (is_retrievable("max_amp", t))
    opt_.max_amp = retrieve_table_double("max_amp", t);
}

void grape::assign_pulse(const sol::table &t) {
  double width = retrieve_table_double("width", t); // unit in ms.
  size_t nsteps = (size_t) (retrieve_table_int("step", t));
  string pattern = "rand_spline";
  if (is_retrievable("init_pattern", t))
    pattern = retrieve_table_str("init_pattern", t);

  double dt = width * 1e-3 / double(nsteps); // into s.

  string str_chs = boost::algorithm::join(superop_.rf_ctrl.chs, " ");
  string code = "user_rf = shapedRF{name = 'opt-rf', width = " + boost::lexical_cast<string>(width) +
      ",  step = " + boost::lexical_cast<string>(nsteps) +
      ",  max_amp = 1" +
      ",  channel = '" + str_chs + "'," +
      "pattern = '" + pattern + "'}";

  g_lua->script(code);

  string s = str(boost::format("pulse width - [%.3f] ms, steps - [%d], step width - [%.3f] us.\n") % width % nsteps
                     % (dt * 1e6));
  ssl_color_text("info", s);

  sol::object val = (*g_lua)["user_rf"];
  seq_block &sb = val.as<seq_block &>(); // note the original type should be 'seq_block'.
  shaped_rf &user_rf = (shaped_rf &) sb; // transfrom into 'shaped_rf'.

  rf_ = &user_rf;
}

void grape::h5write(string file_name) const {
  if (file_name.empty()) {
    string time_s = sys_time();
    file_name = "oc_" + time_s + ".h5";
  }
  H5File file(file_name, H5F_ACC_TRUNC);

  rf_->switch_rf_mode("amp/phase");
  rf_->h5write(file, "opt");

  ssl::utility::h5write(file, nullptr, "obj", stl2vec(opt_.vf));
  file.close();
}

void grape::assign_aux_var() {
  traj_ = state_traj(rf_->get_steps());
  int nbb = superop_.L0s.size();
  if (nbb > 1) {
    omp_set_num_threads(omp_core_num);
    for (size_t i = 0; i < superop_.grad_bb.size(); i++)
      superop_.grad_bb[i] = vector<double>(rf_->get_dims(), 0);

    traj_omp_ = vector<state_traj>(omp_core_num);
    for (int i = 0; i < omp_core_num; i++) {
      traj_omp_[i] = state_traj(rf_->get_steps());
    }
  }
}

void grape::opt_amplitude_constraint(nlopt::opt &opt) {
  if (opt_.max_amp <= 0)
    return;

  size_t dim = rf_->get_dims();
  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  vector<double> up_bound(dim);
  vector<double> low_bound(dim);

  switch (opt_.algo) {
    case nlopt::LD_MMA:
    case nlopt::LD_LBFGS: {
      //if(rf_mode_ ==_amp_phase)
      //      for (size_t i = 0; i < opt_rf_.nsteps; i++) {
      //          for (size_t j = 0; j < opt_rf_.nchannels; j++) {
      //              // amp.
      //              up_bound[2 * opt_rf_.nchannels * i + 2 * j] = opt_rf_.max_amp * 2 * _pi;
      //              low_bound[2 * opt_rf_.nchannels * i + 2 * j] = 0;

      //              // phase.
      //              up_bound[2 * opt_rf_.nchannels * i + 2 * j + 1] = 2 * _pi;
      //              low_bound[2 * opt_rf_.nchannels * i + 2 * j + 1] = 0;
      //          }
      //      }

      if (rf_->mode() == _ux_uy)
        for (size_t i = 0; i < nsteps; i++) {
          for (size_t j = 0; j < nchannels; j++) {
            // ux.
            up_bound[2 * nchannels * i + 2 * j] = opt_.max_amp * 2 * _pi;
            low_bound[2 * nchannels * i + 2 * j] = -opt_.max_amp * 2 * _pi;

            // uy.
            up_bound[2 * nchannels * i + 2 * j + 1] = opt_.max_amp * 2 * _pi;
            low_bound[2 * nchannels * i + 2 * j + 1] = -opt_.max_amp * 2 * _pi;
          }
        }

      opt.set_upper_bounds(up_bound);
      opt.set_lower_bounds(low_bound);
    }
      break;
      /*case nlopt::LD_MMA: {
          amp_constraint_data data = { 0, 0, opt_rf_.nchannels, opt_rf_.max_amp * 2 * _pi };
          for (size_t i = 0; i < opt_rf_.nsteps; i++) {
              for (size_t j = 0; j < opt_rf_.nchannels; j++) {
                  data.i = i;
                  data.j = j;
                  opt.add_inequality_constraint(amplitude_constraint, &data, opt_.opt_f);
              }
          }
      }
      break;*/
    default:break;
  }
}

double grape::objfunc_broadband(const vector<double> &x, vector<double> &grad) {
  rf_->update_raw_data(x.data());
  int N = superop_.L0s.size();
  vec phi = vec::Zero(N);

  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  double dt = rf_->width_in_ms() * 1e-3 / double(nsteps); // into s.

#pragma omp parallel for
  for (int p = 0; p < N; p++) {
    int id = omp_get_thread_num();
    traj_omp_[id].forward[0] = init_state_;
    traj_omp_[id].backward[nsteps] = targ_list_[p];
    //cout << id << "\n";
    sp_cx_mat L;
    sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;
    double kx = 1, ky = 1;
    sp_cx_vec rho = traj_omp_[id].forward[0];
    for (size_t i = 0; i < nsteps; i++) {
      L = L0;
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
      traj_omp_[id].forward[i + 1] = ssl::spinsys::step(rho, L, dt);
      rho = traj_omp_[id].forward[i + 1];
    }
    rho = traj_omp_[id].backward[nsteps];
    for (int i = nsteps - 1; i >= 0; i--) {
      L = L0.adjoint();
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
      traj_omp_[id].backward[i] = ssl::spinsys::step(rho, L, -dt);
      rho = traj_omp_[id].backward[i];
    }
    sp_cx_mat Gx, Gy, tmp;
    int k = 0;
    for (size_t i = 0; i < nsteps; i++) {
      L = L0;
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
      tmp = traj_omp_[id].backward[i + 1].adjoint() * ssl::spinsys::propagator(L, dt);
      for (size_t j = 0; j < nchannels; j++) {
        Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
        Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);

        superop_.grad_bb[p][k] = traced(tmp * Gx * traj_omp_[id].forward[i]).real();
        superop_.grad_bb[p][k + 1] = traced(tmp * Gy * traj_omp_[id].forward[i]).real();
        k += 2;
      }
    }
    phi[p] = transfer_fidelity(traj_omp_[id].forward[nsteps], targ_list_[p]);
    /// transfer_fidelity(targ_list_[p], targ_list_[p]);
  } // end parallel for.

  double val = phi.sum() / (double) N;
  double alpha = 0;

  int k = 0;
  for (size_t i = 0; i < nsteps; i++) {
    for (size_t j = 0; j < nchannels; j++) {
      double gx = 0;
      double gy = 0;
      for (int p = 0; p < N; p++) {
        gx += superop_.grad_bb[p][k];
        gy += superop_.grad_bb[p][k + 1];
      }
      grad[k] = gx / (double) N;
      grad[k + 1] = gy / (double) N;

      //double gxx = grad[k];
      //double gyy = grad[k + 1];

      //if (opt_rf_.rf_mode == _amp_phase)  // Translate the control gradient into phase gradient
      //cartesian2polar(x[k], x[k + 1], gxx, gyy, grad[k], grad[k + 1]);

      // rf power reduction.
      /*double ux, uy;
      ux = x[2 * opt_rf_.nchannels * i + 2 * j];
      uy = x[2 * opt_rf_.nchannels * i + 2 * j + 1];
      grad[k] += 2.0 * alpha*ux*opt_rf_.dt;
      grad[k+1] += 2.0 * alpha*uy*opt_rf_.dt;*/

      k += 2;
    }
  }

  double phi_rf = rf_->average_power();
  double PHI1 = val;
  double PHI2 = alpha * phi_rf * dt;
  //std::cout << boost::format("==> %04d  [%.4f] [%.4f]\n") % (++opt_.iter_count) % PHI1 % PHI2;
  std::cout << boost::format("==> %04d  [%.4f]\n") % (++opt_.iter_count) % PHI1;
  opt_.vf.push_back(PHI1 - PHI2);
  return (PHI1 - PHI2);
}

double grape::objfunc(const vector<double> &x, vector<double> &grad) {
  rf_->update_raw_data(x.data());
  double alpha = 0;
  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  double dt = rf_->width_in_ms() * 1e-3 / double(nsteps); // into s.
  traj_.forward[0] = init_state_;
  traj_.backward[nsteps] = targ_state_;
  sp_cx_mat L;
  sp_cx_mat L0 = superop_.L0 + ci * superop_.R;
  double kx = 1, ky = 1;
  sp_cx_vec rho = traj_.forward[0];
  for (size_t i = 0; i < nsteps; i++) {
    L = L0;
    for (size_t j = 0; j < nchannels; j++)
      L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
    traj_.forward[i + 1] = ssl::spinsys::step(rho, L, dt);
    rho = traj_.forward[i + 1];
  }
  rho = traj_.backward[nsteps];
  for (int i = nsteps - 1; i >= 0; i--) {
    L = L0.adjoint();
    for (size_t j = 0; j < nchannels; j++)
      L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
    traj_.backward[i] = ssl::spinsys::step(rho, L, -dt);
    rho = traj_.backward[i];
  }
  sp_cx_mat Gx, Gy, tmp;
  int k = 0;
  for (size_t i = 0; i < nsteps; i++) {
    L = L0;
    for (size_t j = 0; j < nchannels; j++)
      L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
    tmp = traj_.backward[i + 1].adjoint() * ssl::spinsys::propagator(L, dt);
    for (size_t j = 0; j < nchannels; j++) {
      Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
      Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);
      grad[k] = traced(tmp * Gx * traj_.forward[i]).real();
      grad[k + 1] = traced(tmp * Gy * traj_.forward[i]).real();

      // rf power reduction.
      double ux, uy;
      ux = x[2 * nchannels * i + 2 * j];
      uy = x[2 * nchannels * i + 2 * j + 1];
      grad[k] += 2.0 * alpha * ux * dt;
      grad[k + 1] += 2.0 * alpha * uy * dt;
      k += 2;
    }
  }
  double val = transfer_fidelity(traj_.forward[nsteps], targ_state_);// / transfer_fidelity(targ_state_, targ_state_);
  double phi_rf = rf_->average_power();
  double PHI1 = val;
  double PHI2 = alpha * phi_rf * dt;
  //std::cout << boost::format("==> %04d  [%.4f] [%.4f]\n") % (++opt_.iter_count) % PHI1 % PHI2;
  std::cout << boost::format("==> %04d  [%.4f]\n") % (++opt_.iter_count) % PHI1;
  opt_.vf.push_back(PHI1 - PHI2);
  return (PHI1 - PHI2);
}

void grape::cartesian2polar(double amp, double phase, double gx, double gy, double &g_amp, double &g_phase) {
  double ux = amp * cos(phase);
  double uy = amp * sin(phase);
  g_amp = gx * ux + gy * uy;
  g_amp /= amp;
  g_phase = -gx * uy + gy * ux;
}

double grape::amplitude_constraint(unsigned n, const double *x, double *grad, void *data) {
  amp_constraint_data *d = reinterpret_cast<amp_constraint_data *>(data);

  double ux = x[2 * d->chs * d->i + 2 * d->j];
  double uy = x[2 * d->chs * d->i + 2 * d->j + 1];

  return ux * ux + uy * uy - d->amp * d->amp;
}

sp_cx_mat grape::update_rf_ham(const double *x,
                               size_t step,
                               size_t channel,
                               size_t nchannels,
                               double kx,
                               double ky) const {
  double ux = x[2 * nchannels * step + 2 * channel];
  double uy = x[2 * nchannels * step + 2 * channel + 1];

  /*if (opt_rf_.rf_mode ==_amp_phase) {
      double uxx = ux*cos(uy);
      double uyy = ux*sin(uy);
      ux = uxx;
      uy = uyy;
  }*/
  // Rf inhom
  ux *= kx;
  uy *= ky;
  return ux * superop_.rf_ctrl.Lx[channel] + uy * superop_.rf_ctrl.Ly[channel];
}
sp_cx_mat grape::propagator_derivative(const sp_cx_mat &L, const sp_cx_mat &Lrf, double dt) {
  sp_cx_mat commu1, commu2, commu3;
  commu1 = commutator(L, Lrf);
  commu2 = commutator(L, commu1);
  commu3 = commutator(L, commu2);
  return -ci * dt * Lrf + 1 / 2.0 * dt * dt * commu1 + 1 / 6.0 * dt * dt * dt * commu2
      - 1 / 24.0 * dt * dt * dt * dt * commu3;
}

void grape::projection(const sol::table &t) {
  sol::object val = retrieve_table("init_state", t);
  sp_cx_vec init_state = sys_->smart_state(val.as<string>());
  init_state = norm_state(init_state);
  //init_state = levante_ernst(init_state);

  val = retrieve_table("rf", t);
  seq_block &sb = val.as<seq_block &>();
  shaped_rf &user_rf = (shaped_rf &) sb; // transfrom into 'shaped_rf'.
  user_rf.convert2(_ux_uy);

  size_t nsteps = user_rf.get_steps();
  size_t nchannels = user_rf.get_channels();
  double dt = user_rf.width_in_ms() * 1e-3 / double(nsteps); // into s.

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
    //obsrv_state.push_back(levante_ernst(iter->second));
    expr.push_back(iter->first);
  }

  string s;
  for (size_t p = 0; p < expr.size(); p++)
    s += "row " + to_string(p + 1) + "-" + expr[p] + " ";
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

  vector<double> x = user_rf.clone_raw_data();
  if (str_opt == "step") {
    sp_cx_mat L;
    sp_cx_mat L0 = superop_.L0 + superop_.R;
    sp_cx_vec *forward = new sp_cx_vec[nsteps + 1];
    forward[0] = init_state;
    sp_cx_vec rho = forward[0];

    for (size_t i = 0; i < nsteps; i++) {
      L = L0;
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(x.data(), i, j, nchannels);
      forward[i + 1] = ssl::spinsys::step(rho, L, dt);
      rho = forward[i + 1];

      for (size_t k = 0; k < obsrv_state.size(); k++) {
        sp_cx_vec compo = obsrv_state[k];
        transfer_wave(i, k) = transfer_fidelity(rho, compo);// / transfer_fidelity(compo, compo);
      }
    }

    // TO BE REMOVED (for NSFC project).
    /*
    for (size_t i = 0; i < nsteps / 2; i++) {
      L = L0;
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(x.data(), i, j, nchannels);
      forward[i + 1] = ssl::spinsys::step(rho, L, dt);
      rho = forward[i + 1];

      for (size_t k = 0; k < obsrv_state.size(); k++) {
        sp_cx_vec compo = obsrv_state[k];
        transfer_wave(i, k) = transfer_fidelity(rho, compo);// / transfer_fidelity(compo, compo);
      }
    }

    //cout<< (superop_.rf_ctrl.Lx[0]*superop_.L0-superop_.L0*superop_.rf_ctrl.Lx[0])<<"\n";
    //exit(0);
    sp_cx_mat Lrf = superop_.rf_ctrl.Lx[0] + superop_.rf_ctrl.Lx[1];
    sp_cx_vec tmp = ssl::spinsys::step(rho, Lrf, _pi);
    rho = tmp;

    for (size_t i = nsteps / 2; i < nsteps; i++) {
      L = L0;
      for (size_t j = 0; j < nchannels; j++)
        L += update_rf_ham(x.data(), i, j, nchannels);
      forward[i + 1] = ssl::spinsys::step(rho, L, dt);
      rho = forward[i + 1];

      for (size_t k = 0; k < obsrv_state.size(); k++) {
        sp_cx_vec compo = obsrv_state[k];
        transfer_wave(i, k) = transfer_fidelity(rho, compo);// / transfer_fidelity(compo, compo);
      }
    }
     */
    // TO BE REMOVED (for NSFC project).

  } else if (str_opt == "broadband") {
    omp_set_num_threads(omp_core_num);

    vector<sp_cx_vec *> forward_trajs_parfor = vector<sp_cx_vec *>(omp_core_num);
    for (int i = 0; i < omp_core_num; i++) {
      sp_cx_vec *forward = new sp_cx_vec[nsteps + 1];
      forward[0] = init_state;
      forward_trajs_parfor[i] = forward;
    }

#pragma omp parallel for
    for (int p = 0; p < nrows; p++) {
      int id = omp_get_thread_num();

      sp_cx_mat L;
      sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;
      sp_cx_vec rho = forward_trajs_parfor[id][0];

      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, nchannels);
        forward_trajs_parfor[id][i + 1] = ssl::spinsys::step(rho, L, dt);
        rho = forward_trajs_parfor[id][i + 1];
      }

      // TO BE REMOVED (for NSFC project).
      /*
      for (size_t i = 0; i < nsteps / 2; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, nchannels);
        forward_trajs_parfor[id][i + 1] = ssl::spinsys::step(rho, L, dt);
        rho = forward_trajs_parfor[id][i + 1];
      }

      sp_cx_mat Lrf = superop_.rf_ctrl.Lx[0] + superop_.rf_ctrl.Lx[1];
      sp_cx_vec tmp = ssl::spinsys::step(rho, Lrf, _pi);
      rho = tmp;

      for (size_t i = nsteps / 2; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, nchannels);
        forward_trajs_parfor[id][i + 1] = ssl::spinsys::step(rho, L, dt);
        rho = forward_trajs_parfor[id][i + 1];
      }
       */
      // TO BE REMOVED (for NSFC project).

      for (size_t k = 0; k < obsrv_state.size(); k++) {
        sp_cx_vec compo = obsrv_state[k];
        transfer_wave(p, k) = transfer_fidelity(rho, compo);// / transfer_fidelity(compo, compo);
      }
    }
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
    time = vec::LinSpaced(nrows, 0, user_rf.width_in_ms());
    fig_spec = "title<transfer trajectories of basis operators> xlabel<pulse duration / ms> ylabel<coefficient>";
  } else if (str_opt == "broadband") {
    int n = superop_.nominal_offset.size();
    time = vec::LinSpaced(nrows, superop_.nominal_offset[0], superop_.nominal_offset[n - 1]);
    //fig_spec = "title[transfer coefficient] xlabel[freq offset / kHz]";
    fig_spec = "title<transfer coefficient> xlabel<J coupling / Hz>";
  }

  if (expr.size() > 5)
    fig_spec += " gnuplot<set key outside>";

  fig_spec += " color<Spectral>";
  //fig_spec += " latex[ON]";

  string lege;
  for (size_t i = 0; i < expr.size(); i++)
    lege += expr[i] + ";";
  fig_spec += "legend<" + lege + ">";

  plot(fig_spec, line_series(time, lines));

  ofstream ofstr("basis_projection.dat");
  ofstr.precision(3);
  ofstr << transfer_wave;
  ofstr.close();

  return; // Currently ignore the following map plot.

  vec2 xrange, yrange;
  if (str_opt == "step") {
    xrange[0] = 0;
    xrange[1] = user_rf.width_in_ms();
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
  //g_lua->script("plot('title#transfer coefficient# xlabel#pulse duration / ms# ylabel#projection component# gnuplot#set label front \"ttt\" at graph 0.5,0.5#', _map)");
  if (str_opt == "step")
    g_lua->script(
        "plot('title<transfer coefficient> xlabel<pulse duration / ms> ylabel<projection component> gnuplot<set label front \""
            + s + "\"  at graph 0.5,0.5 center font \",10\">', _map)");
  if (str_opt == "broadband")
    g_lua->script(
        "plot('title<transfer coefficient> xlabel<J coupling / Hz> ylabel<projection component> gnuplot<set label front \""
            + s + "\"  at graph 0.5,0.5 center font \",10\">', _map)");
}

} /* namespace oc */
} /* namespace ssl */
