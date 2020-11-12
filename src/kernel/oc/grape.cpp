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
#include <chrono>
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
  sys_ = &sys;
  superop_.rf_ctrl = sys.rf_hamiltonian();
  superop_.L0 = sys.free_hamiltonian();
  superop_.R = sys.relaxation();
  superop_.L0s = sys.free_hamiltonians();
  int bb_size = superop_.L0s.size();
  //std::cout << bb_size << "\n";
  if (bb_size) {
    // superop_.profile = vec::Ones(bb_size);
    superop_.nominal_offset = sys.nominal_broadband();
    superop_.grad_bb = std::vector<std::vector<double>>(bb_size);
  } else {  // no offset case.
    superop_.L0s.push_back(superop_.L0);
    superop_.grad_bb = std::vector<std::vector<double>>(1);
  }
  opt_model_ = _rho2rho;
  axis_ = uxuy_;
}

grape::~grape() {
}

double grape::maxf() const {
  return max_val_;
}

seq_block &grape::optimize(const sol::table &t) {
  assign_pulse(t);
  assign_x();
  assign_nlopt(t);
  assign_state(t);
  assign_aux_var();

  if (opt_model_ == _rho2rho) 
      optimizer_->set_max_objective(objfunc_broadband, this);
  if (opt_model_ == _propagator)
    optimizer_->set_max_objective(objfunc_propagator, this);

  print();
  double max_f;
  nlopt::result result = optimizer_->optimize(x_, max_f);
  if (result == nlopt::SUCCESS)
    ssl_color_text("info", "pulse optimization succeed.\n");
  if (result == nlopt::FAILURE)
    ssl_color_text("err", "pulse optimization failed.\n");
  if (result == nlopt::MAXEVAL_REACHED)
    ssl_color_text("warn", "pulse optimization terminated due to maximum iterations limitation.\n");

  max_val_ = max_f;
  rf_->update_raw_data(axis_, x_.data());

  grape::h5write();
  return *rf_;
}

double grape::objfunc_broadband(const std::vector<double> &x, std::vector<double> &grad,
                                void *func) {
  return ((grape *)func)->objfunc_broadband(x, grad);
}
double grape::objfunc_propagator(const std::vector<double> &x, std::vector<double> &grad,
                                 void *func) {
  return ((grape *)func)->objfunc_propagator(x, grad);
}

void grape::assign_state(const sol::table &t) {
  if (is_retrievable("targ_op", t)) {
    opt_model_ = _propagator;
    targ_op_ = sys_->smart_op(retrieve_table_str("targ_op", t));
    targ_op_ = ssl::spinsys::propagator(targ_op_, _pi);
    int n = targ_op_.rows();
    init_op_ = sp_cx_mat(n, n);
    init_op_.setIdentity();

    // in case of unified targ for diff freq offsets.
    if (superop_.L0s.size())
      targ_op_list_ = std::vector<sp_cx_mat>(superop_.L0s.size(), targ_op_);
    return;
  }

  opt_model_ = _rho2rho;
  init_state_ = sys_->smart_state(retrieve_table_str("init_state", t));
  init_state_ = norm_state(init_state_);

  sol::object obj = retrieve_table("targ_state", t);
  if (obj.get_type() == sol::type::string) {
    targ_state_ = sys_->smart_state(retrieve_table_str("targ_state", t));
    targ_state_ = norm_state(targ_state_);

    // in case of unified targ for diff freq offsets.
    if (superop_.L0s.size())
      targ_list_ = std::vector<sp_cx_vec>(superop_.L0s.size(), targ_state_);
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
  // TO BE REMOVE.
  epsilon_ = 0;      // default val.
  epsilon_num_ = 1;  // default val.
  if (is_retrievable("epsilon", t)) {
    epsilon_ = retrieve_table_double("epsilon", t);
  }

  if (is_retrievable("epsilon_num", t)) {
    epsilon_num_ = retrieve_table_int("epsilon_num", t);
  }
  nlopt::algorithm alg = nlopt::LD_LBFGS;

  if (is_retrievable("algorithm", t)) {
    std::string oc = retrieve_table_str("algorithm", t);
    boost::to_upper(oc);
    if (oc == "MNA") alg = nlopt::LD_MMA;
  }
  optimizer_ = new nlopt::opt(alg, x_dim_);

  if (is_retrievable("max_eval", t)) {
    int val = retrieve_table_int("max_eval", t);
    if (val > 0)
      optimizer_->set_maxeval(val);
    else
      throw std::runtime_error("max_eval must > 0.");
  }

  if (is_retrievable("xtol_rel", t)) {
    double val = retrieve_table_double("xtol_rel", t);
    optimizer_->set_xtol_rel(val);
  } else
    optimizer_->set_xtol_rel(1e-8);

  assign_constraint(t);
}
void grape::print() const {
  if (optimizer_->get_algorithm() == nlopt::LD_MMA)
    ssl_color_text("oc", "nlopt algorithm used: LD_MMA.\n");
  if (optimizer_->get_algorithm() == nlopt::LD_LBFGS)
    ssl_color_text("oc", "nlopt algorithm used: LD_LBFGS.\n");
  
  ssl_color_text("oc", "x dim: " + boost::lexical_cast<std::string>(x_dim_) +  ".\n");
  ssl_color_text("oc", "nlopt xtol_rel: " + boost::lexical_cast<std::string>(optimizer_->get_xtol_rel()) +  ".\n");

  
}
void grape::assign_x() {
  x_.clear();
  rf_->convert2(_ux_uy);  // must do it.

  if (axis_ == uxuy_)
    x_ = rf_->clone_raw_data();
  else if (axis_ == ux_)
    x_ = rf_->clone_raw_data_ux();
  else if (axis_ == uy_)
    x_ = rf_->clone_raw_data_uy();

  x_dim_ = x_.size();
}

void grape::assign_pulse(const sol::table &t) {
  double width = retrieve_table_double("width", t); // unit in ms.
  size_t nsteps = (size_t) (retrieve_table_int("step", t));
  std::string pattern = "rand_spline";
  if (is_retrievable("init_pattern", t))
    pattern = retrieve_table_str("init_pattern", t);

  double max_init_amp = 1;
  if (is_retrievable("max_init_amp", t))
    max_init_amp = retrieve_table_double("max_init_amp", t);


  double dt = width * 1e-3 / double(nsteps); // into s.

  std::string str_chs;
  if (is_retrievable("limit_channel", t))
      str_chs = retrieve_table_str("limit_channel", t);
  else
  str_chs = boost::algorithm::join(superop_.rf_ctrl.chs, " ");

  std::string code = "user_rf = shapedRF{name = 'opt-rf', width = " +
                boost::lexical_cast<std::string>(width) +
                ",  step = " + boost::lexical_cast<std::string>(nsteps) +
                ",  max_amp = " + boost::lexical_cast<std::string>(max_init_amp) +
                ",  channel = '" + str_chs + "'," + "pattern = '" + pattern +
                "'}";

  g_lua->script(code);

   std::string s = str(boost::format("pulse width - [%.3f] ms, steps - [%d], step width - [%.3f] us.\n") % width % nsteps
                     % (dt * 1e6));
  ssl_color_text("info", s);

  sol::object val = (*g_lua)["user_rf"];
  seq_block &sb = val.as<seq_block &>(); // note the original type should be 'seq_block'.
  shaped_rf &user_rf = (shaped_rf &) sb; // transfrom into 'shaped_rf'.

  rf_ = &user_rf;

  if (is_retrievable("limit_axis", t)) {
    std::string s = retrieve_table_str("limit_axis", t);
    if (s == "x")
      axis_ = ux_;
    else if (s == "y")
      axis_ = uy_;
  }

  rf_->convert2(_ux_uy);  // must do it.

  if(axis_==ux_) {
   std::vector<double>zeros = std::vector<double>(rf_->get_dims() / 2, 0);
    rf_->update_raw_data_uy(zeros.data());
  }
  if (axis_ == uy_) {
    std::vector<double> zeros = std::vector<double>(rf_->get_dims() / 2, 0);
    rf_->update_raw_data_ux(zeros.data());
  }
}

void grape::h5write(std::string file_name) const {
  if (file_name.empty()) {
    std::string time_s = sys_time();
    file_name = "oc_" + time_s + ".h5";
  }
  H5File file(file_name, H5F_ACC_TRUNC);

  rf_->switch_rf_mode("amp/phase");
  rf_->h5write(file, "opt");

  ssl::utility::h5write(file, nullptr, "obj", stl2vec(obj_val_));
  file.close();
}

void grape::assign_aux_var() {
  // traj_ = state_traj(rf_->get_steps());
  int nbb = superop_.L0s.size();
  if (nbb >= 1) {
    omp_set_num_threads(omp_core_num);
    for (size_t i = 0; i < superop_.grad_bb.size(); i++)
      superop_.grad_bb[i] = std::vector<double>(rf_->get_dims(), 0);

    if (opt_model_ == _rho2rho) {
      traj_omp_ = std::vector<state_traj>(omp_core_num);
      for (int i = 0; i < omp_core_num; i++) {
        traj_omp_[i] = state_traj(rf_->get_steps());
      }
    }
    if (opt_model_ == _propagator) {
      op_traj_omp_ = std::vector<op_traj>(omp_core_num);
      for (int i = 0; i < omp_core_num; i++) {
        op_traj_omp_[i] = op_traj(rf_->get_steps());
      }
    }
  }
}

void grape::assign_constraint(const sol::table &t) {
  if (is_retrievable("max_amp", t)) {
    double val = retrieve_table_double("max_amp", t);  // unit in hz.
    val *= 2 * _pi;

    size_t dim = rf_->get_dims();
    size_t nsteps = rf_->get_steps();
    size_t nchannels = rf_->get_channels();
    std::vector<double> up_bound(dim);
    std::vector<double> low_bound(dim);

    // only for ux/uy mode.
    if (rf_->mode() == _ux_uy)
      for (size_t i = 0; i < nsteps; i++) {
        for (size_t j = 0; j < nchannels; j++) {
          // ux.
          up_bound[2 * nchannels * i + 2 * j] = val;
          low_bound[2 * nchannels * i + 2 * j] = -val;

          // uy.
          up_bound[2 * nchannels * i + 2 * j + 1] = val;
          low_bound[2 * nchannels * i + 2 * j + 1] = -val;
        }
      }

    optimizer_->set_upper_bounds(up_bound);
    optimizer_->set_lower_bounds(low_bound);
  }
}

double grape::objfunc_broadband(const std::vector<double> &x, std::vector<double> &grad) {
  auto start = std::chrono::system_clock::now();
  rf_->update_raw_data(axis_, x.data());
  int N = superop_.L0s.size();
  vec phi = vec::Zero(N);

  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  double dt = rf_->get_dt() * 1e-6;  // into s.
  std::vector<std::string> chs = rf_->get_channels_str();

  // double alpha = alpha0;
  vec rf_scaling =
      vec::LinSpaced(epsilon_num_, -epsilon_, epsilon_);  // e.g. -10% - 10%

#pragma omp parallel for
  for (int p = 0; p < N; p++) {
    int id = omp_get_thread_num();
    traj_omp_[id].forward[0] = init_state_;
    traj_omp_[id].backward[nsteps] = targ_list_[p];
    sp_cx_mat L;
    sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;

    int dim = rf_->get_dims();
    if (axis_ != uxuy_) dim /= 2;
    superop_.grad_bb[p] = std::vector<double>(dim, 0);
    // start rf inhom.

    for (int q = 0; q < rf_scaling.size();
         q++) {  // for each rf scaling factor, do

      double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];

      sp_cx_vec rho = traj_omp_[id].forward[0];
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
        traj_omp_[id].forward[i + 1] = ssl::spinsys::step(rho, L, dt);
        rho = traj_omp_[id].forward[i + 1];
      }
      rho = traj_omp_[id].backward[nsteps];
      for (int i = nsteps - 1; i >= 0; i--) {
        L = L0.adjoint();
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
        traj_omp_[id].backward[i] = ssl::spinsys::step(rho, L, -dt);
        rho = traj_omp_[id].backward[i];
      }
      sp_cx_mat Gx, Gy, tmp;
      int k = 0;
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
        tmp = traj_omp_[id].backward[i + 1].adjoint() *
              ssl::spinsys::propagator(L, dt);
        for (size_t j = 0; j < nchannels; j++) {
          if (axis_ == uxuy_) {
            Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
            Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);

            superop_.grad_bb[p][k] +=
                traced(tmp * Gx * traj_omp_[id].forward[i]).real();
            superop_.grad_bb[p][k + 1] +=
                traced(tmp * Gy * traj_omp_[id].forward[i]).real();
            k += 2;
          } else if (axis_ == ux_) {
            Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);

            superop_.grad_bb[p][k] +=
                traced(tmp * Gx * traj_omp_[id].forward[i]).real();
            k += 1;
          } else if (axis_ == uy_) {
            Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);

            superop_.grad_bb[p][k] +=
                traced(tmp * Gy * traj_omp_[id].forward[i]).real();
            k += 1;
          }
        }
      }
      phi[p] += transfer_fidelity(traj_omp_[id].forward[nsteps], targ_list_[p]);
    }

    /// transfer_fidelity(targ_list_[p], targ_list_[p]);
  }  // end parallel for.

  double val = phi.sum() / (double)N / (double)rf_scaling.size();

  int k = 0;
  for (size_t i = 0; i < nsteps; i++) {
    for (size_t j = 0; j < nchannels; j++) {
      double gx = 0;
      double gy = 0;
      if (axis_ == uxuy_) {
        for (int p = 0; p < N; p++) {
          gx += superop_.grad_bb[p][k];
          gy += superop_.grad_bb[p][k + 1];
        }
        grad[k] = gx / (double)N / (double)rf_scaling.size();
        grad[k + 1] = gy / (double)N / (double)rf_scaling.size();

        // rf power reduction.
        /*if (alpha != 0) {
          double ux, uy;
          ux = x[2 * nchannels * i + 2 * j];
          uy = x[2 * nchannels * i + 2 * j + 1];
          grad[k] -= 2.0 * alpha * ux * dt;
          grad[k + 1] -= 2.0 * alpha * uy * dt;
        }*/
        k += 2;
      } else if (axis_ == ux_) {
        for (int p = 0; p < N; p++) {
          gx += superop_.grad_bb[p][k];
        }
        grad[k] = gx / (double)N / (double)rf_scaling.size();

        k += 1;
      } else if (axis_ == uy_) {
        for (int p = 0; p < N; p++) {
          gy += superop_.grad_bb[p][k];
        }
        grad[k] = gy / (double)N / (double)rf_scaling.size();

        k += 1;
      }
    }
  }
  auto end = std::chrono::system_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  //double PHI1 = val;
  //double PHI2 = 0;
  // double PHI2 = alpha * rf_->rf_power();
  // std::cout << boost::format("==> %04d  [%.8f] [%.8f]\n") %
  // (++opt_.iter_count) % PHI1 % PHI2;
  std::cout << boost::format("==> %04d  [%.8f]\n") % (optimizer_->get_numevals()) % val;
 obj_val_.push_back(val);
  //std::cout << "Use Time:" << double(duration.count()) * microseconds::period::num / microseconds::period::den << " s.\n";
  return val;
}

double grape::objfunc_propagator(const std::vector<double> &x,
                                 std::vector<double> &grad) {
  rf_->update_raw_data(x.data());
  int N = superop_.L0s.size();
  vec phi = vec::Zero(N);

  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  double dt = rf_->get_dt() * 1e-6;  // into s.
  std::vector<std::string> chs = rf_->get_channels_str();

  // double alpha = alpha0;
  vec rf_scaling =
      vec::LinSpaced(epsilon_num_, -epsilon_, epsilon_);  // e.g. -10% - 10%

#pragma omp parallel for
  for (int p = 0; p < N; p++) {
    int id = omp_get_thread_num();
    op_traj_omp_[id].forward[0] = init_op_;
    op_traj_omp_[id].backward[nsteps] = targ_op_list_[p];
    sp_cx_mat L;
    sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;

     int dim = rf_->get_dims();
    if (axis_ != uxuy_) dim /= 2;
    superop_.grad_bb[p] = std::vector<double>(dim, 0);
    // start rf inhom.

    for (int q = 0; q < rf_scaling.size();
         q++) {  // for each rf scaling factor, do

      double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];

      sp_cx_mat rho = op_traj_omp_[id].forward[0];
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
        op_traj_omp_[id].forward[i + 1] = ssl::spinsys::propagator(L, dt) * rho;
        rho = op_traj_omp_[id].forward[i + 1];
      }
      rho = op_traj_omp_[id].backward[nsteps];
      for (int i = nsteps - 1; i >= 0; i--) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
        op_traj_omp_[id].backward[i] =
            ssl::spinsys::propagator(L, dt).adjoint() * rho;
        rho = op_traj_omp_[id].backward[i];
      }
      sp_cx_mat Gx, Gy, tmp;
      int k = 0;
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
        tmp = op_traj_omp_[id].backward[i + 1].adjoint() *
              ssl::spinsys::propagator(L, dt);
        for (size_t j = 0; j < nchannels; j++) {
          if (axis_ == uxuy_) {
            Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
            Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);

            superop_.grad_bb[p][k] +=
                traced(tmp * Gx * op_traj_omp_[id].forward[i]).real();
            superop_.grad_bb[p][k + 1] +=
                traced(tmp * Gy * op_traj_omp_[id].forward[i]).real();
            k += 2;
          } else if (axis_ == ux_) {
            Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);

            superop_.grad_bb[p][k] +=
                traced(tmp * Gx * op_traj_omp_[id].forward[i]).real();
            k += 1;
          } else if (axis_ == uy_) {
            Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);
            superop_.grad_bb[p][k] +=
                traced(tmp * Gy * op_traj_omp_[id].forward[i]).real();
            k += 2;
          }
        }
      }
      phi[p] += transfer_fidelity(op_traj_omp_[id].forward[nsteps],
                                  targ_op_list_[p]) /
                transfer_fidelity(targ_op_list_[p], targ_op_list_[p]);
    }

    /// transfer_fidelity(targ_list_[p], targ_list_[p]);
  }  // end parallel for.

  double val = phi.sum() / (double)N / (double)rf_scaling.size();

  int k = 0;
  for (size_t i = 0; i < nsteps; i++) {
    for (size_t j = 0; j < nchannels; j++) {
      double gx = 0;
      double gy = 0;
     if (axis_ == uxuy_) {
        for (int p = 0; p < N; p++) {
          gx += superop_.grad_bb[p][k];
          gy += superop_.grad_bb[p][k + 1];
        }
        grad[k] = gx / (double)N / (double)rf_scaling.size();
        grad[k + 1] = gy / (double)N / (double)rf_scaling.size();

        // rf power reduction.
        /*if (alpha != 0) {
          double ux, uy;
          ux = x[2 * nchannels * i + 2 * j];
          uy = x[2 * nchannels * i + 2 * j + 1];
          grad[k] -= 2.0 * alpha * ux * dt;
          grad[k + 1] -= 2.0 * alpha * uy * dt;
        }*/
        k += 2;
      } else if (axis_ == ux_) {
        for (int p = 0; p < N; p++) {
          gx += superop_.grad_bb[p][k];
        }
        grad[k] = gx / (double)N / (double)rf_scaling.size();

        k += 1;
      } else if (axis_ == uy_) {
        for (int p = 0; p < N; p++) {
          gy += superop_.grad_bb[p][k];
        }
        grad[k] = gy / (double)N / (double)rf_scaling.size();

        k += 1;
      }
    }
  }

  //double PHI1 = val;
  //double PHI2 = 0;
  // double PHI2 = alpha * rf_->rf_power();
  // std::cout << boost::format("==> %04d  [%.8f] [%.8f]\n") %
  // (++opt_.iter_count) % PHI1 % PHI2;
  std::cout << boost::format("==> %04d  [%.8f]\n") % (optimizer_->get_numevals()) % val;
  obj_val_.push_back(val);
  return (val);
}
// void grape::cartesian2polar(double amp, double phase, double gx, double gy,
//                            double &g_amp, double &g_phase) {
//  double ux = amp * cos(phase);
//  double uy = amp * sin(phase);
//  g_amp = gx * ux + gy * uy;
//  g_amp /= amp;
//  g_phase = -gx * uy + gy * ux;
//}

double grape::amplitude_constraint(unsigned n, const double *x, double *grad,
                                   void *data) {
  amp_constraint_data *d = reinterpret_cast<amp_constraint_data *>(data);

  double ux = x[2 * d->chs * d->i + 2 * d->j];
  double uy = x[2 * d->chs * d->i + 2 * d->j + 1];

  return ux * ux + uy * uy - d->amp * d->amp;
}

sp_cx_mat grape::update_rf_ham(const double *x, size_t step, size_t channel, std::string ch_str,
                               size_t nchannels, double kx, double ky) {
  int ch = superop_.rf_ctrl.channel_index(ch_str);
  if (axis_ == uxuy_) {
    double ux = x[2 * nchannels * step + 2 * channel];
    double uy = x[2 * nchannels * step + 2 * channel + 1];
    // Rf inhom
    ux *= kx;
    uy *= ky;
    return ux * superop_.rf_ctrl.Lx[ch] + uy * superop_.rf_ctrl.Ly[ch];
  } else if (axis_ == ux_) {
    double ux = x[nchannels * step + channel];
    // Rf inhom
    ux *= kx;
    return ux * superop_.rf_ctrl.Lx[ch];
  } else if (axis_ == uy_) {
    double uy = x[nchannels * step + channel];
    // Rf inhom
    uy *= ky;
    return uy * superop_.rf_ctrl.Ly[ch];
  }
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
  epsilon_ = 0;
  epsilon_num_ = 1;  // default val.
  if (is_retrievable("epsilon", t)) {
    epsilon_ = retrieve_table_double("epsilon", t);
  }
  if (is_retrievable("epsilon_num", t)) {
    epsilon_num_ = retrieve_table_int("epsilon_num", t);
  }
  sol::object val = retrieve_table("init_state", t);
  sp_cx_vec init_state = sys_->smart_state(val.as<std::string>());
  init_state = norm_state(init_state);
  //init_state = levante_ernst(init_state);

  val = retrieve_table("rf", t);
  seq_block &sb = val.as<seq_block &>();
 shaped_rf &user_rf = (shaped_rf &) sb; // transfrom into 'shaped_rf'.
  user_rf.convert2(_ux_uy);

  size_t nsteps = user_rf.get_steps();
  size_t nchannels = user_rf.get_channels();
  double dt = user_rf.width_in_ms() * 1e-3 / double(nsteps); // into s.

  std::map<std::string, sp_cx_vec> obsrv_state_map;
  std::vector<std::string> expr;
  std::vector<sp_cx_vec> obsrv_state;
  if (is_retrievable("observ_states", t)) {
    val = retrieve_table("observ_states", t);
    sol::table expr_table = val.as<sol::table>();
    if (expr_table.size() == 0) {
      obsrv_state_map = sys_->cartesian_basis_states();
    } else {
      for (size_t i = 0; i < expr_table.size(); i++) {
        sol::object val = expr_table[i + 1];
        std::string exp = val.as<std::string>();
        sp_cx_vec rho = sys_->smart_state(exp);
        obsrv_state_map.insert(std::pair<std::string, sp_cx_vec>(exp, rho));
      }
    }

  } else if (is_retrievable("ignore_states", t)) {
    val = retrieve_table("ignore_states", t);
    sol::table expr_table = val.as<sol::table>();
    if (expr_table.size() == 0) {
      obsrv_state_map = sys_->cartesian_basis_states();
    } else {
      obsrv_state_map = sys_->cartesian_basis_states();

      for (size_t i = 0; i < expr_table.size(); i++) {
        sol::object val = expr_table[i + 1];
        std::string exp = val.as<std::string>();
        std::map<std::string, sp_cx_vec>::iterator key = obsrv_state_map.find(exp);
        if (key != obsrv_state_map.end()) obsrv_state_map.erase(key);
      }
    }
  } else {
    obsrv_state_map = sys_->cartesian_basis_states();
  }

  std::map<std::string, sp_cx_vec>::iterator iter;
  for (iter = obsrv_state_map.begin(); iter != obsrv_state_map.end(); iter++) {
    obsrv_state.push_back(norm_state(iter->second));
    // obsrv_state.push_back(levante_ernst(iter->second));
    expr.push_back(iter->first);
  }
  std::string dim1 = "observed states\n";
  std::string dim2 = "";
  std::string dim3 = "rf scaling\n";
  std::string s;
  for (size_t p = 0; p < expr.size(); p++)
    s += std::to_string(p + 1) + "-" + expr[p] + "\n";
  dim1 += s;

  std::string str_opt = "step";
  if (is_retrievable("option", t)) str_opt = retrieve_table_str("option", t);

  if (str_opt != "step" && str_opt != "broadband") {
    std::string s = "unknown projection option ** " + str_opt +
               " ** using 'step' or 'broadband' instead.";
    throw std::runtime_error(s.c_str());
  }

  std::vector<double> x = user_rf.clone_raw_data();
  vec rf_scaling =
      vec::LinSpaced(epsilon_num_, -epsilon_, epsilon_);  // e.g. -10% - 10%

  size_t n = rf_scaling.size();
  if (n == 1)
    dim3 += "no rf inhomogeneity";
  else
    dim3 += boost::lexical_cast<std::string>(rf_scaling[0] * 100) + ":" +
            boost::lexical_cast<std::string>(rf_scaling[n - 1] * 100) + " " +
            boost::lexical_cast<std::string>(n) + " %\n";

  cube comp_dist;
  // for BROADBAND case: states, freq offsets, rf scalings
  // for STEP case: states, steps, rf scalings
  std::vector<std::string> chs = rf_->get_channels_str();

  if (str_opt == "step") {
    dim2 = "pulse steps\n";
    dim2 += "interval: " + boost::lexical_cast<std::string>(dt) + " s\n";
    comp_dist = cube(obsrv_state.size(), nsteps,
                     rf_scaling.size());  // states, steps, rf scalings
    sp_cx_mat L;
    sp_cx_mat L0 = superop_.L0 + ci * superop_.R;
    sp_cx_vec *forward = new sp_cx_vec[nsteps + 1];
    forward[0] = init_state;

    for (int q = 0; q < rf_scaling.size();
         q++) {  // for each rf scaling factor, do
      double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];
      sp_cx_vec rho = forward[0];

      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
        forward[i + 1] = ssl::spinsys::step(rho, L, dt);
        rho = forward[i + 1];

        for (size_t k = 0; k < obsrv_state.size(); k++) {
          sp_cx_vec compo = obsrv_state[k];
          comp_dist(k, i, q) = transfer_fidelity(rho, compo);
        }
      }
    }
  } else if (str_opt == "broadband") {
    dim2 = "freq offsets\n";
    int n = superop_.nominal_offset.size();
    dim2 += boost::lexical_cast<std::string>(superop_.nominal_offset[0]) + ":" +
            boost::lexical_cast<std::string>(superop_.nominal_offset[n - 1]) + " " +
            boost::lexical_cast<std::string>(n) + " Hz\n";

    comp_dist = cube(obsrv_state.size(), superop_.L0s.size(),
                     rf_scaling.size());  // states, freq offsets, rf scalings

    omp_set_num_threads(omp_core_num);

    std::vector<sp_cx_vec *> forward_trajs_parfor =
        std::vector<sp_cx_vec *>(omp_core_num);
    for (int i = 0; i < omp_core_num; i++) {
      sp_cx_vec *forward = new sp_cx_vec[nsteps + 1];
      forward[0] = init_state;
      forward_trajs_parfor[i] = forward;
    }

#pragma omp parallel for
    for (int p = 0; p < superop_.L0s.size(); p++) {
      int id = omp_get_thread_num();

      for (int q = 0; q < rf_scaling.size();
           q++) {  // for each rf scaling factor, do
        double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];
        sp_cx_mat L;
        sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;
        sp_cx_vec rho = forward_trajs_parfor[id][0];

        for (size_t i = 0; i < nsteps; i++) {
          L = L0;
          for (size_t j = 0; j < nchannels; j++)
            L += update_rf_ham(x.data(), i, j, chs[j], nchannels, kx, ky);
          forward_trajs_parfor[id][i + 1] = ssl::spinsys::step(rho, L, dt);
          rho = forward_trajs_parfor[id][i + 1];
        }

        for (size_t k = 0; k < obsrv_state.size(); k++) {
          sp_cx_vec compo = obsrv_state[k];
          // transfer_wave(p, k) = transfer_fidelity(rho, compo);  // /
          // transfer_fidelity(compo, compo);
          comp_dist(k, p, q) = transfer_fidelity(rho, compo);
        }
      }
    }
  }

  std::string time_s = sys_time();
  H5File file("proj_" + time_s + ".h5", H5F_ACC_TRUNC);
  ssl::utility::h5write(file, nullptr, "projection", comp_dist);
  ssl::utility::h5write(file, nullptr, "dim1", dim1);
  ssl::utility::h5write(file, nullptr, "dim2", dim2);
  ssl::utility::h5write(file, nullptr, "dim3", dim3);
  file.close();

  if (rf_scaling.size() == 1) {  // no rf inhomogeneity.
    Eigen::Tensor<double, 2> sub =
        comp_dist.chip(0, 2);  // rf scaling at dim 2.
    Eigen::Map<mat> m(sub.data(), sub.dimension(0), sub.dimension(1));
    mat grid = m.matrix();

    sol::table lines;

    lines = g_lua->create_table();
    for (int i = 0; i < grid.rows(); i++) {  // each state.
      vec line = grid.row(i).transpose();
      lines.add(line);
    }

    std::string fig_spec;
    vec xval;
    if (str_opt == "step") {
      xval = vec::LinSpaced(grid.cols(), 0, user_rf.width_in_ms()); // transfer trajectories of basis operators
	  std::string s = retrieve_table_str("init_state", t);
      fig_spec = "title<initial state: " + s + "> ";
      fig_spec +=
          "xlabel<pulse " //initial state I_{1y}
          "duration "
          "/ ms> ylabel<magnetization>";
      fig_spec +=
          "xrange<0:" + boost::lexical_cast<std::string>(user_rf.width_in_ms()) +
          "> ";
    } else if (str_opt == "broadband") {
      int n = superop_.nominal_offset.size();
      xval = vec::LinSpaced(grid.cols(), superop_.nominal_offset[0],
                            superop_.nominal_offset[n - 1]);
      xval *= 1e-3;
      fig_spec = "xlabel<frequency offset / kHz> ylabel<magnetization>"; // title<scan 1> 
      fig_spec += "xrange<" + boost::lexical_cast<std::string>(xval[0]) + ":" +
                  boost::lexical_cast<std::string>(xval[xval.size() - 1]) + "> ";
    }
    if (expr.size() > 5)
        fig_spec += " gnuplot<set key outside>";
      //fig_spec += " gnuplot<set label 'Scan 2' at graph 0.5,0.5 center font 'Arial,26'\n set ytics 0.2\n set key horizontal above>";
	  //fig_spec += " gnuplot<set label 'Scan 2' at graph 0.5,0.5 center font 'Arial,26'\n set ytics 0.2\n unset key>";
	  //fig_spec += " gnuplot<set ytics 0.2\n unset key>";

    fig_spec += " lw<7>";
    fig_spec += " color<YiZhang16,16>";
    std::string lege;
    for (size_t i = 0; i < expr.size(); i++) lege += expr[i] + ";";
    fig_spec += " legend<" + lege + ">";
    plot(fig_spec, line_series(xval, lines));
    return;
  }

  // rf scaling case.

  Eigen::Tensor<double, 2> sub =
      comp_dist.chip(0, 0);  // observed states at dim 0, by default, plot the
                             // map for first state.
  std::string state_name = expr[0];
  Eigen::Map<mat> m(sub.data(), sub.dimension(0), sub.dimension(1));
  mat grid = m.matrix();

  vec2 xrange, yrange;
  yrange[0] = superop_.nominal_offset[0];
  yrange[1] = superop_.nominal_offset[superop_.nominal_offset.size() - 1];
  yrange *= 1e-3;  // into kHz.
  xrange[0] = rf_scaling[0];
  xrange[1] = rf_scaling[rf_scaling.size() - 1];
  xrange *= 100;  // into %

  utility::map gnu(grid, "style<3d>");
  gnu.xrange = xrange;
  gnu.yrange = yrange;
  std::string fig_spec =
      "'ylabel<freq offset / kHz> xlabel<rf inhomogeneity / %> color<Spectral> "
      "gnuplot<set xlabel offset -1,-1; set ylabel offset -2,-2; set palette "
      "negative; set zrange [0.8:1]; set cbrange [0.8:1]> ";

  (*g_lua)["_comp_dist"] = gnu;
  g_lua->script("plot(" + fig_spec + " title<magnetization map [" + state_name +
                "]>', _comp_dist)");
}

} /* namespace oc */
} /* namespace ssl */
