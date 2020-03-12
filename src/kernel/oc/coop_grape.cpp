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
#include <chrono>
using namespace chrono;
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;

namespace ssl {
namespace oc {
coop_model coop_grape::coop_model_;  // by default, use multi-scan cooperative
                                     // pulses optimization model.

struct ms_coop_ctrl {
  sp_cx_vec rho_bar;
  sp_cx_vec rho_lambda;
  vector<sp_cx_vec> obsrv_state;
  vec beta;
  vec gamma;
  void alloc(spin_system *sys, sp_cx_vec targ_state) {
    std::map<string, sp_cx_vec> obsrv_state_map;
    obsrv_state_map = sys->cartesian_basis_states();
    std::map<string, sp_cx_vec>::iterator iter;
    for (iter = obsrv_state_map.begin(); iter != obsrv_state_map.end();
         iter++) {
      obsrv_state.push_back(norm_state(iter->second));
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

  void calc_rho_bar(const vector<state_traj>& ms_traj, size_t nsteps/*, sp_cx_vec targ_state*/) {
	   /*for (size_t i = 0; i < obsrv_state.size(); i++)
      gamma[i] = transfer_fidelity(targ_state, obsrv_state[i]);*/

    rho_bar.setZero();

    for (int i = 0; i < ms_traj.size(); i++)
      rho_bar += ms_traj[i].forward[nsteps];

    rho_bar /= (double)ms_traj.size(); // average.
    rho_bar.coeffRef(0) = 1;  // IMPROTANT! DO NOT NORMALIED FOR BB CASE.

    //rho_bar.coeffRef(0) = 0;
    //rho_bar = norm_state(rho_bar); // need to be normalizd FOR NORMAL CASE.
  }
  void calc_lambdaN(size_t ncoop) {
    for (size_t i = 0; i < obsrv_state.size(); i++)
      beta[i] = transfer_fidelity(rho_bar, obsrv_state[i]);

    rho_lambda.setZero();
    for (size_t i = 0; i < obsrv_state.size(); i++)
      rho_lambda += (gamma[i] - beta[i]) * obsrv_state[i];
    rho_lambda *= 2.0/ double(ncoop);
    rho_lambda.coeffRef(0) = 1;  // NOTE.
    // IF rho_lambda NEED TO BE normalizd??
    //rho_lambda.coeffRef(0) = 0;
    //rho_lambda = norm_state(rho_lambda);
  }

  double calc_obj() {
    vec val = gamma - beta;
    return 1- val.cwiseAbs2().sum();
  }
};

coop_grape::coop_grape(spin_system &sys) : grape(sys) { is_coop_ = true; coop_model_=_ms_coop;}

coop_grape::~coop_grape() {}

void coop_grape::assign_pulse(const sol::table &t) {
  grape::assign_pulse(t);
  size_t ncoop = retrieve_table_size_t("ncoop", t);
  for (size_t i = 0; i < ncoop; i++) {
    shaped_rf *ptr = rf_->Clone();
    ptr->copy_config_table();
    ptr->set_name("ms-coop " + to_string(i + 1) + "-" + to_string(ncoop));
    ptr->assign();
    coop_rf_.push_back(ptr);
  }

  coop_model_ = _ms_coop;
  if (is_retrievable("coop_model", t)) {
    string s = retrieve_table_str("coop_model", t);
    if (s == "ss-coop") coop_model_ = _ss_coop;
  }

  if (is_retrievable("insert_op", t)) {
    int sys_dim = superop_.L0.rows();
    sp_cx_mat Ie = sp_cx_mat(sys_dim, sys_dim);
    Ie.setIdentity();
    inserted_ops_ =
        vector<sp_cx_mat>(coop_rf_.size() - 1, Ie);  // note size ncoop-1.

    sol::object obj = retrieve_table("insert_op", t);
    sol::table par_table = obj.as<sol::table>();

    sol::object val = par_table[1];  // 1st element: operator.
    sp_cx_mat op = sys_->smart_op(val.as<string>());
    op = ssl::spinsys::propagator(op, _pi);
    val = par_table[2];  // 2nd element: id(1,2,3...)
    inserted_ops_[val.as<size_t>() - 1] = op;

    coop_model_ = _ss_coop;
  }
}

void coop_grape::assign_aux_var() {
  int nbb = superop_.L0s.size();
  coop_par_ = new ms_coop_ctrl;
  coop_par_->alloc(sys_, targ_list_[0]);

  for (size_t i = 0; i < coop_rf_.size(); i++)
    coop_traj_.push_back(state_traj(rf_->get_steps()));

  omp_set_num_threads(omp_core_num);
}

sol::object coop_grape::optimize(const sol::table &t, sol::this_state s) {
  assign_state(t);  // shared with grape.
  assign_nlopt(t);  // shared with grape.
  assign_pulse(t);
  assign_aux_var();

  vector<double> x;
  for (size_t i = 0; i < coop_rf_.size(); i++) {
    coop_rf_[i]->convert2(_ux_uy);
    vector<double> raw = coop_rf_[i]->clone_raw_data();
    // ONLY FOR TEST.
    vector<double> raw_x;
    for (int i=0;i<raw.size();i++)
      if (i % 2 == 0) raw_x.push_back(raw[i]);

    if (i==1) {
      for (size_t j = 0; j < raw.size(); j++) {
        raw[i] *= -1;
      }
    }
    x.insert(x.end(), raw_x.begin(), raw_x.end());  // column major.
    raw_x.clear();
    //x.insert(x.end(), raw.begin(), raw.end());  // column major.
    raw.clear();
  }

  // optimization processing.
  nlopt::opt opt(opt_.algo, x.size());
  opt.set_xtol_rel(opt_.tol_f);
  if (opt_.max_eval > 0) opt.set_maxeval(opt_.max_eval);
  if (opt_.max_time > 0) opt.set_maxtime(opt_.max_time);
  if (opt_.stopval > 0) opt.set_stopval(opt_.stopval);

  if (opt_model_ == _rho2rho) 
      opt.set_max_objective(objfunc_broadband, this);

  // boundary limitation.
  opt_amplitude_constraint(opt);

  double max_f = 0;
  nlopt::result result = opt.optimize(x, max_f);
  if (result == nlopt::SUCCESS)
    ssl_color_text("info", "pulse optimization succeed.\n");
  if (result == nlopt::FAILURE)
    ssl_color_text("err", "pulse optimization failed.\n");
  if (result == nlopt::MAXEVAL_REACHED)
    ssl_color_text("warn",
                   "pulse optimization terminated due to maximum iterations "
                   "limitation.\n");

  opt_.max_f = max_f;


  // ONLY FOR TEST.
  int dimx = (int)rf_->get_dims()/2;
  vector<double> xx;
  for (int i=0;i<dimx;i++)
  {
    xx.push_back(x[i]);
    xx.push_back(0);
  }
  for (int i = 0; i < dimx; i++) {
    xx.push_back(x[i+dimx]);
     xx.push_back(0);
  }
  cout << xx.size() << "\n";
  sol::state_view lua(s);
  sol::table rfs = lua.create_table();
  int dim = (int)rf_->get_dims();
  for (size_t i = 0; i < coop_rf_.size(); i++) {
    coop_rf_[i]->update_raw_data(xx.data() + i * dim);
    rfs.add((seq_block *)coop_rf_[i]);
  }
  h5write();
  return rfs;
}

void coop_grape::opt_amplitude_constraint(nlopt::opt &opt) {
  if (opt_.max_amp <= 0) return;

  switch (opt_.algo) {
    case nlopt::LD_MMA:
    case nlopt::LD_LBFGS: {

		size_t ncoop = coop_rf_.size();
	  size_t dim = rf_->get_dims();
      size_t nsteps = rf_->get_steps();
      size_t nchannels = rf_->get_channels();
      vector<double> up_bound(dim * ncoop);
      vector<double> low_bound(dim * ncoop);

      for (size_t k = 0; k < coop_rf_.size(); k++) {
        if (coop_rf_[k]->mode() == _ux_uy)
          for (size_t i = 0; i < nsteps; i++) {
            for (size_t j = 0; j < nchannels; j++) {
              // ux.
              up_bound[dim*k + 2 * nchannels * i + 2 * j] = opt_.max_amp * 2 * _pi;
              low_bound[dim*k + 2 * nchannels * i + 2 * j] = -opt_.max_amp * 2 * _pi;

              // uy.
              up_bound[dim*k + 2 * nchannels * i + 2 * j + 1] = opt_.max_amp * 2 * _pi;
              low_bound[dim * k + 2 * nchannels * i + 2 * j + 1] = -opt_.max_amp * 2 * _pi;
            }
          }
      }
      opt.set_upper_bounds(up_bound);
      opt.set_lower_bounds(low_bound);
    } break;
    default:
      break;
  }
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

double coop_grape::objfunc_broadband(const vector<double> &x, vector<double> &g,
                                     void *func) {
  if (coop_model_ == _ms_coop)
   return ((coop_grape *)func)->objfunc_broadband_ms_coop(x, g);  // co_objfunc
  else
  return ((coop_grape *)func)->objfunc_broadband_ss_coop(x, g); 
}
double coop_grape::objfunc_broadband_ss_coop(const vector<double> &x,
                                             vector<double> &g) {
  auto start = system_clock::now();
  int N = superop_.L0s.size();
  vec phi = vec::Zero(N);

  int ncoop = (int)coop_rf_.size();

  vector<Eigen::Map<const vec>> shapes;  // pointer.
  vector<Eigen::Map<vec>> grads;         // pointer.
  vector<vec> grads_ss_coop;
  int idx = 0;
  for (int i = 0; i < ncoop; i++) {
    int dim = (int)coop_rf_[i]->get_dims()/2; // ONLY FOR TEST. /2 ONLY X AXIS.
    Eigen::Map<const vec> shape(x.data() + idx, dim);
    Eigen::Map<vec> grad(g.data() + idx, dim);
    grad.setZero();
    shapes.push_back(shape);
    grads.push_back(grad);
    vec g(dim);
    g.setZero();
    grads_ss_coop.push_back(g);
    idx += dim;
  }

  vector<vector<vec>> grad_bb(N, grads_ss_coop);

#pragma omp parallel for
  for (int p = 0; p < N; p++) {
    sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;

    vector<state_traj> co_traj_rho;
    for (int index = 0; index < ncoop; index++) {
      state_traj cur_index = state_traj(coop_rf_[index]->get_steps());

      if (index == 0)  // only for 1st pulse since it is ss-coop.
        cur_index.forward[0] = init_state_;
      if (index == ncoop - 1)  // only for last pulse since it is ss-coop.
        cur_index.backward[coop_rf_[index]->get_steps()] = targ_list_[p];

      co_traj_rho.push_back(cur_index);
    }
    //cout << "##1\n";
    //#pragma omp parallel for
    for (int index = 0; index < ncoop; index++) {
      sp_cx_mat L;
      double kx = 1, ky = 1;
      size_t nsteps = coop_rf_[index]->get_steps();
      size_t nchannels = coop_rf_[index]->get_channels();
      double dt =
          coop_rf_[index]->width_in_ms() * 1e-3 / double(nsteps);  // into s.

      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(shapes[index], i, j, nchannels, kx, ky);
        co_traj_rho[index].forward[i + 1] =
            ssl::spinsys::step(co_traj_rho[index].forward[i], L, dt);
      }

      // init for the next pulse.
      if (index + 1 < ncoop)  // do not do for the last pulse.
        co_traj_rho[index + 1].forward[0] =
            inserted_ops_[index] * co_traj_rho[index].forward[nsteps];
    }
    //cout << "##2\n";
    //#pragma omp parallel for
    for (int index = ncoop - 1; index >= 0; index--) {
      sp_cx_mat L;
      double kx = 1, ky = 1;
      size_t nsteps = coop_rf_[index]->get_steps();
      size_t nchannels = coop_rf_[index]->get_channels();
      double dt =
          coop_rf_[index]->width_in_ms() * 1e-3 / double(nsteps);  // into s.

      for (int i = nsteps - 1; i >= 0; i--) {
        L = L0.adjoint();
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(shapes[index], i, j, nchannels, kx, ky);
        co_traj_rho[index].backward[i] =
            ssl::spinsys::step(co_traj_rho[index].backward[i + 1], L, -dt);
      }

      // init for the previous pulse.
      if (index - 1 >= 0)  // do not do for the first pulse.
        co_traj_rho[index - 1].backward[coop_rf_[index - 1]->get_steps()] =
            inserted_ops_[index - 1].adjoint() * co_traj_rho[index].backward[0];
    }
    //cout << "##3\n";
    //#pragma omp parallel for
    for (int index = 0; index < ncoop; index++) {
      sp_cx_mat Gx, Gy, tmp;
      sp_cx_mat L;
      int k = 0;
      double kx = 1, ky = 1;
      size_t nsteps = coop_rf_[index]->get_steps();
      size_t nchannels = coop_rf_[index]->get_channels();
      double dt =
          coop_rf_[index]->width_in_ms() * 1e-3 / double(nsteps);  // into s.
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(shapes[index], i, j, nchannels, kx, ky);
        tmp = co_traj_rho[index].backward[i + 1].adjoint() *
              ssl::spinsys::propagator(L, dt);
        for (size_t j = 0; j < nchannels; j++) {
          Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
          //Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);
          grad_bb[p][index](k) =
              traced(tmp * Gx * co_traj_rho[index].forward[i]).real();
          //grad_bb[p][index](k + 1) =
          //    traced(tmp * Gy * co_traj_rho[index].forward[i]).real();
          //k += 2; // ONLY FOR TEST.
          k += 1;
        }
      }
    }
    //cout << "##4\n";
    phi[p] = transfer_fidelity(
        co_traj_rho[ncoop - 1].forward[coop_rf_[ncoop - 1]->get_steps()],
        targ_list_[p]);

	for (size_t i = 0; i < co_traj_rho.size(); i++)
		co_traj_rho[i].free();
  }

  for (int index = 0; index < ncoop; index++) {
    for (int p = 0; p < grad_bb.size(); p++) grads[index] += grad_bb[p][index];
    grads[index] /= (double)N;
  }

  double val = phi.sum() / (double)N;

  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  std::cout << boost::format("==> %04d [%.7f]\n") % (++opt_.iter_count) % val;
  //cout << "Use Time:" << double(duration.count()) * microseconds::period::num / microseconds::period::den << " s.\n";
  opt_.vf.push_back(val);
  return val;
}

double coop_grape::objfunc_broadband_ms_coop(const vector<double> &x,
                                             vector<double> &g) {
  int ncoop = (int)coop_rf_.size();
  int dim = (int)rf_->get_dims();
  Eigen::Map<const mat> shape(x.data(), dim, ncoop);
  Eigen::Map<mat> grad(g.data(), dim, ncoop);
  grad.setZero();

  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  double dt = rf_->width_in_ms() * 1e-3 / double(nsteps);  // into s.

  int N = superop_.L0s.size();
  vec phi = vec::Zero(N);

  mat mg = mat(dim, ncoop);
  mg.setZero();
  vector<mat> grad_bb(N, mg);

#pragma omp parallel for
  for (int p = 0; p < N; p++) {
    sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;

    vector<state_traj> co_traj_rho;
    for (int index = 0; index < ncoop; index++) {
      state_traj cur_index = state_traj(nsteps);
      cur_index.forward[0] = init_state_;
      co_traj_rho.push_back(cur_index);
    }

    //#pragma omp parallel for
    for (int index = 0; index < ncoop; index++) {
      sp_cx_mat L;
      double kx = 1, ky = 1;
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(shape, index, i, j, nchannels, kx, ky);
        co_traj_rho[index].forward[i + 1] =
            ssl::spinsys::step(co_traj_rho[index].forward[i], L, dt);
      }
    }

    ms_coop_ctrl my = *coop_par_;
    my.calc_rho_bar(co_traj_rho, nsteps/*, targ_list_[p]*/);
    my.calc_lambdaN(ncoop);

    phi[p] = my.calc_obj();

    //#pragma omp parallel for
    for (int index = 0; index < ncoop; index++) {
      co_traj_rho[index].backward[nsteps] = my.rho_lambda;

      sp_cx_mat L;
      double kx = 1, ky = 1;
      for (int i = nsteps - 1; i >= 0; i--) {
        L = L0.adjoint();
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(shape, index, i, j, nchannels, kx, ky);
        co_traj_rho[index].backward[i] =
            ssl::spinsys::step(co_traj_rho[index].backward[i + 1], L, -dt);
      }
    }

    //#pragma omp parallel for
    for (int index = 0; index < ncoop; index++) {
      sp_cx_mat Gx, Gy, tmp;
      sp_cx_mat L;
      int k = 0;
      double kx = 1, ky = 1;
      for (size_t i = 0; i < nsteps; i++) {
        L = L0;
        for (size_t j = 0; j < nchannels; j++)
          L += update_rf_ham(shape, index, i, j, nchannels, kx, ky);
        tmp = co_traj_rho[index].backward[i + 1].adjoint() *
              ssl::spinsys::propagator(L, dt);
        for (size_t j = 0; j < nchannels; j++) {
          Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
          Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);
          grad_bb[p](k, index) =
              traced(tmp * Gx * co_traj_rho[index].forward[i]).real();
          grad_bb[p](k + 1, index) =
              traced(tmp * Gy * co_traj_rho[index].forward[i]).real();
          k += 2;
        }
      }
    }

		for (size_t i = 0; i < co_traj_rho.size(); i++)
		co_traj_rho[i].free();
  }

  for (int i = 0; i < grad_bb.size(); i++) grad += grad_bb[i];

  grad /= (double)N;
  double val = phi.sum() / (double)N;

  std::cout << boost::format("==> %04d [%.7f]\n") % (++opt_.iter_count) % val;
  opt_.vf.push_back(val);
  return val;
}

double coop_grape::co_objfunc(const vector<double> &x, vector<double> &g) {
  int ncoop = (int)coop_rf_.size();
  int dim = (int)rf_->get_dims();
  Eigen::Map<const mat> shape(x.data(), dim, ncoop);
  Eigen::Map<mat> grad(g.data(), dim, ncoop);

  size_t nsteps = rf_->get_steps();
  size_t nchannels = rf_->get_channels();
  double dt = rf_->width_in_ms() * 1e-3 / double(nsteps);  // into s.

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

  coop_par_->calc_rho_bar(coop_traj_, nsteps/*, targ_list_[0]*/);
  coop_par_->calc_lambdaN(ncoop);

#pragma omp parallel for
  for (int p = 0; p < ncoop; p++) {
    coop_traj_[p].backward[nsteps] = coop_par_->rho_lambda;

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
      tmp = coop_traj_[p].backward[i + 1].adjoint() *
            ssl::spinsys::propagator(L, dt);
      for (size_t j = 0; j < nchannels; j++) {
        Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
        Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);
        grad(k, p) = traced(tmp * Gx * coop_traj_[p].forward[i]).real();
        grad(k + 1, p) = traced(tmp * Gy * coop_traj_[p].forward[i]).real();
        k += 2;
      }
    }
  }

  double val = coop_par_->calc_obj();
  std::cout << boost::format("==> %04d [%.7f]\n") % (++opt_.iter_count) % val;
  opt_.vf.push_back(val);
  return val;
}

void coop_grape::projection(const sol::table &t) {
  epsilon_ = 0;
  epsilon_num_ = 1;  // default val.
  if (is_retrievable("epsilon", t)) {
    epsilon_ = retrieve_table_double("epsilon", t);
  }
  if (is_retrievable("epsilon_num", t)) {
    epsilon_num_ = retrieve_table_int("epsilon_num", t);
  }
  sol::object val = retrieve_table("init_state", t);
  sp_cx_vec init_state = sys_->smart_state(val.as<string>());
  init_state = norm_state(init_state);
  // init_state = levante_ernst(init_state);

  val = retrieve_table("coop_rf", t);
  sol::table rf_table = val.as<sol::table>();

  vector<shaped_rf *> co_rfs;
  for (size_t i = 0; i < rf_table.size(); i++) {
    sol::object val = rf_table[i + 1];
    seq_block *p_rf = val.as<seq_block *>();
    co_rfs.push_back((shaped_rf *)p_rf);
    // co_rfs[i]->switch_rf_mode("ux/uy");
  }

  size_t nsteps = co_rfs[0]->get_steps();
  size_t nchannels = co_rfs[0]->get_channels();
  double dt = co_rfs[0]->width_in_ms() * 1e-3 / double(nsteps);  // into s.

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

  } else if (is_retrievable("ignore_states", t)) {
    val = retrieve_table("ignore_states", t);
    sol::table expr_table = val.as<sol::table>();
    if (expr_table.size() == 0) {
      obsrv_state_map = sys_->cartesian_basis_states();
    } else {
      obsrv_state_map = sys_->cartesian_basis_states();

      for (size_t i = 0; i < expr_table.size(); i++) {
        sol::object val = expr_table[i + 1];
        string exp = val.as<string>();
        std::map<string, sp_cx_vec>::iterator key = obsrv_state_map.find(exp);
        if (key != obsrv_state_map.end()) obsrv_state_map.erase(key);
      }
    }
  } else {
    obsrv_state_map = sys_->cartesian_basis_states();
  }

  std::map<string, sp_cx_vec>::iterator iter;
  for (iter = obsrv_state_map.begin(); iter != obsrv_state_map.end(); iter++) {
    obsrv_state.push_back(norm_state(iter->second));
    // obsrv_state.push_back(levante_ernst(iter->second));
    expr.push_back(iter->first);
  }
  string dim1 = "observed states\n";
  string dim2 = "";
  string dim3 = "rf scaling\n";
  string s;
  for (size_t p = 0; p < expr.size(); p++)
    s += to_string(p + 1) + "-" + expr[p] + "\n";
  dim1 += s;

  string str_opt = "step";
  if (is_retrievable("option", t)) str_opt = retrieve_table_str("option", t);

  if (str_opt != "step" && str_opt != "broadband") {
    string s = "unknown projection option ** " + str_opt +
               " ** using 'step' or 'broadband' instead.";
    throw std::runtime_error(s.c_str());
  }

  vec rf_scaling =
      vec::LinSpaced(epsilon_num_, -epsilon_, epsilon_);  // e.g. -10% - 10%

  size_t n = rf_scaling.size();
  if (n == 1)
    dim3 += "no rf inhomogeneity";
  else
    dim3 += boost::lexical_cast<string>(rf_scaling[0] * 100) + ":" +
            boost::lexical_cast<string>(rf_scaling[n - 1] * 100) + " " +
            boost::lexical_cast<string>(n) + " %\n";

  cube comp_dist;
  // for BROADBAND case: states, freq offsets, rf scalings
  // for STEP case: states, steps, rf scalings
  //  if given 'insert_op', the coop mode must be ss-coop.
	    if (is_retrievable("insert_op", t)) {
      inserted_ops_.clear();
      int sys_dim = superop_.L0.rows();
      sp_cx_mat Ie = sp_cx_mat(sys_dim, sys_dim);
      Ie.setIdentity();
      inserted_ops_ =
          vector<sp_cx_mat>(co_rfs.size() - 1, Ie);  // note size ncoop-1.

      sol::object obj = retrieve_table("insert_op", t);
      sol::table par_table = obj.as<sol::table>();

      sol::object val = par_table[1];  // 1st element: operator.
      sp_cx_mat op = sys_->smart_op(val.as<string>());
      op = ssl::spinsys::propagator(op, _pi);
      val = par_table[2];  // 2nd element: id(1,2,3...)
      inserted_ops_[val.as<size_t>() - 1] = op;

      coop_model_ = _ss_coop;
            }


	if (coop_model_ == _ss_coop) {
    vector<double> x;
    for (size_t i = 0; i < co_rfs.size(); i++) {
      co_rfs[i]->convert2(_ux_uy);
      vector<double> raw = co_rfs[i]->clone_raw_data();
      x.insert(x.end(), raw.begin(), raw.end());  // column major.
      raw.clear();
    }

    int ncoop = (int)co_rfs.size();

    vector<Eigen::Map<const vec>> shapes;  // pointer.
    int idx = 0;
    for (int i = 0; i < ncoop; i++) {
      int dim = (int)co_rfs[i]->get_dims();
      Eigen::Map<const vec> shape(x.data() + idx, dim);
      shapes.push_back(shape);
      idx += dim;
    }
    if (str_opt == "step") {
      dim2 = "pulse steps\n";
      dim2 += "interval: " + boost::lexical_cast<string>(dt) + " s\n";

      int ss_nsteps = 0;
      for (int index = 0; index < ncoop; index++)
        ss_nsteps += co_rfs[index]->get_steps();

      comp_dist = cube(obsrv_state.size(), ss_nsteps,
                       rf_scaling.size());  // states, steps, rf scalings

	  cout << ss_nsteps << "\n";
      sp_cx_mat L0 = superop_.L0 + ci * superop_.R;

      for (int q = 0; q < rf_scaling.size();
           q++) {  // for each rf scaling factor, do
        double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];

        vector<state_traj> co_traj_rho;
        for (size_t i = 0; i < co_rfs.size(); i++)
          co_traj_rho.push_back(state_traj(co_rfs[i]->get_steps()));

        // init for the 1st pulse.
        co_traj_rho[0].forward[0] = init_state;

        int ss_nsteps = 0;

        for (int index = 0; index < ncoop; index++) {
          sp_cx_mat L;
          size_t nsteps = co_rfs[index]->get_steps();
          size_t nchannels = co_rfs[index]->get_channels();
          double dt = co_rfs[index]->width_in_ms() * 1e-3 /
                      double(nsteps);  // into s.

          for (size_t i = 0; i < nsteps; i++) {
            L = L0;
            for (size_t j = 0; j < nchannels; j++)
              L += update_rf_ham(shapes[index], i, j, nchannels, kx, ky);
            co_traj_rho[index].forward[i + 1] =
                ssl::spinsys::step(co_traj_rho[index].forward[i], L, dt);

            for (size_t k = 0; k < obsrv_state.size(); k++) {
              sp_cx_vec compo = obsrv_state[k];
              comp_dist(k, ss_nsteps + (int)i, q) =
                  transfer_fidelity(co_traj_rho[index].forward[i + 1], compo);
            }

          }
          ss_nsteps += nsteps;
          // init for the next pulse.
          if (index + 1 < ncoop)  // do not do for the last pulse.
            co_traj_rho[index + 1].forward[0] =
                inserted_ops_[index] * co_traj_rho[index].forward[nsteps];
        }
      }
    } else if (str_opt == "broadband") {
      dim2 = "freq offsets\n";
      int n = superop_.nominal_offset.size();
      dim2 += boost::lexical_cast<string>(superop_.nominal_offset[0]) + ":" +
              boost::lexical_cast<string>(superop_.nominal_offset[n - 1]) +
              " " + boost::lexical_cast<string>(n) + " Hz\n";

      comp_dist = cube(obsrv_state.size(), superop_.L0s.size(),
                       rf_scaling.size());  // states, freq offsets, rf scalings
      comp_dist.setZero();
      omp_set_num_threads(omp_core_num);
#pragma omp parallel for
      for (int p = 0; p < superop_.L0s.size(); p++) {
        for (int q = 0; q < rf_scaling.size();
             q++) {  // for each rf scaling factor, do
          double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];

          vector<state_traj> co_traj_rho;
          for (size_t i = 0; i < co_rfs.size(); i++)
            co_traj_rho.push_back(state_traj(co_rfs[i]->get_steps()));

          // init for the 1st pulse.
          co_traj_rho[0].forward[0] = init_state;

          int ss_nsteps = 0;

          sp_cx_mat L;
          sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;

          for (int index = 0; index < ncoop; index++) {
            sp_cx_mat L;
            size_t nsteps = co_rfs[index]->get_steps();
            size_t nchannels = co_rfs[index]->get_channels();
            double dt = co_rfs[index]->width_in_ms() * 1e-3 /
                        double(nsteps);  // into s.

            for (size_t i = 0; i < nsteps; i++) {
              L = L0;
              for (size_t j = 0; j < nchannels; j++)
                L += update_rf_ham(shapes[index], i, j, nchannels, kx, ky);
              co_traj_rho[index].forward[i + 1] =
                  ssl::spinsys::step(co_traj_rho[index].forward[i], L, dt);
            }

            // init for the next pulse.
            if (index + 1 < ncoop)  // do not do for the last pulse.
              co_traj_rho[index + 1].forward[0] =
                  inserted_ops_[index] * co_traj_rho[index].forward[nsteps];
          }

          for (size_t k = 0; k < obsrv_state.size(); k++) {
            sp_cx_vec compo = obsrv_state[k];
            comp_dist(k, p, q) = transfer_fidelity(
                co_traj_rho[ncoop - 1].forward[co_rfs[ncoop - 1]->get_steps()],
                compo);
          }
        }
      }
    }
  } else
  if (coop_model_ == _ms_coop) {
    vector<double> x;
    for (size_t i = 0; i < co_rfs.size(); i++) {
      co_rfs[i]->convert2(_ux_uy);
      vector<double> raw = co_rfs[i]->clone_raw_data();
      x.insert(x.end(), raw.begin(), raw.end());  // column major.
      raw.clear();
    }
    int ncoop = co_rfs.size();
    Eigen::Map<const mat> shape(x.data(), co_rfs[0]->get_dims(), ncoop);

    if (str_opt == "step") {
      dim2 = "pulse steps\n";
      dim2 += "interval: " + boost::lexical_cast<string>(dt) + " s\n";
      comp_dist = cube(obsrv_state.size(), nsteps,
                       rf_scaling.size());  // states, steps, rf scalings
      sp_cx_mat L0 = superop_.L0 + ci * superop_.R;

      for (int q = 0; q < rf_scaling.size();
           q++) {  // for each rf scaling factor, do
        double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];

        vector<state_traj> co_traj_rho;
        for (size_t i = 0; i < co_rfs.size(); i++)
          co_traj_rho.push_back(state_traj(nsteps));

        for (int p = 0; p < ncoop; p++) co_traj_rho[p].forward[0] = init_state;

        for (size_t i = 0; i < nsteps; i++) {
          sp_cx_mat L;

          for (int p = 0; p < ncoop; p++) {
            L = L0;

            sp_cx_vec rho = co_traj_rho[p].forward[i];

            for (size_t j = 0; j < nchannels; j++)
              L += update_rf_ham(shape, p, i, j, nchannels, kx, ky);

            co_traj_rho[p].forward[i + 1] = ssl::spinsys::step(rho, L, dt);
          }

          sp_cx_vec rho_bar = init_state;
          rho_bar.setZero();

          for (int p = 0; p < ncoop; p++)
            rho_bar += co_traj_rho[p].forward[i + 1];

          rho_bar.coeffRef(0) = 0;
          rho_bar = norm_state(rho_bar);

          for (size_t k = 0; k < obsrv_state.size(); k++) {
            sp_cx_vec compo = obsrv_state[k];
            comp_dist(k, i, q) = transfer_fidelity(rho_bar, compo);
          }
        }
      }
    } else if (str_opt == "broadband") {
      dim2 = "freq offsets\n";
      int n = superop_.nominal_offset.size();
      dim2 += boost::lexical_cast<string>(superop_.nominal_offset[0]) + ":" +
              boost::lexical_cast<string>(superop_.nominal_offset[n - 1]) +
              " " + boost::lexical_cast<string>(n) + " Hz\n";

      comp_dist = cube(obsrv_state.size(), superop_.L0s.size(),
                       rf_scaling.size());  // states, freq offsets, rf scalings
      comp_dist.setZero();
      omp_set_num_threads(omp_core_num);
#pragma omp parallel for
      for (int p = 0; p < superop_.L0s.size(); p++) {
        for (int q = 0; q < rf_scaling.size();
             q++) {  // for each rf scaling factor, do
          double kx = 1 + rf_scaling[q], ky = 1 + rf_scaling[q];

          vector<state_traj> co_traj_rho;
          for (size_t i = 0; i < co_rfs.size(); i++)
            co_traj_rho.push_back(state_traj(nsteps));

          for (int index = 0; index < ncoop; index++)
            co_traj_rho[index].forward[0] = init_state;

          sp_cx_mat L;
          sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;

          for (size_t i = 0; i < nsteps; i++) {
            for (int index = 0; index < ncoop; index++) {
              sp_cx_mat L = L0;
              sp_cx_vec rho = co_traj_rho[index].forward[i];
              for (size_t j = 0; j < nchannels; j++)
                L += update_rf_ham(shape, index, i, j, nchannels, kx, ky);

            co_traj_rho[index].forward[i + 1] = ssl::spinsys::step(rho, L, dt);
            }
          }

          sp_cx_vec rho_bar = init_state;
          rho_bar.setZero();
          for (int index = 0; index < ncoop; index++)
            rho_bar += co_traj_rho[index].forward[nsteps];
          rho_bar *= 1 / (double)ncoop;
          rho_bar.coeffRef(0) = 1;

          for (size_t k = 0; k < obsrv_state.size(); k++) {
            sp_cx_vec compo = obsrv_state[k];
            comp_dist(k, p, q) = transfer_fidelity(rho_bar, compo);
          }
        }
      }
    }
  }
  string time_s = sys_time();
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

    string fig_spec;
    vec xval;
    if (str_opt == "step") {
      double ms = 0;
      if (coop_model_ == _ms_coop) ms = co_rfs[0]->width_in_ms();
      if (coop_model_ == _ss_coop) { 
        for (int index = 0; index < co_rfs.size(); index++)
          ms += co_rfs[index]->width_in_ms();  // may be not suit for variable
                                               // dt of different pulses.
      }
      cout << ms << " " << grid.cols() << "\n";
      xval = vec::LinSpaced(grid.cols(), 0, ms);  // transfer trajectories of basis operators
      fig_spec =
          "title<initial state I_{1x}+I_{1y}> xlabel<pulse "
          "duration "
          "/ ms> ylabel<magnetization>";
      fig_spec +=
          "xrange<0:" + boost::lexical_cast<string>(ms) +
          "> ";
    } else if (str_opt == "broadband") {
      int n = superop_.nominal_offset.size();
      xval = vec::LinSpaced(grid.cols(), superop_.nominal_offset[0],
                            superop_.nominal_offset[n - 1]);
      fig_spec = "xlabel<J / Hz> ylabel<magnetization>";
	  //xval *= 1e-3;
      //fig_spec = "xlabel<frequency offset / kHz> ylabel<magnetization>";  // transfer
                                                                   // coefficient
                                                                   // title<2-scan
                                                                   // MS-COOP>
      fig_spec += "xrange<" + boost::lexical_cast<string>(xval[0]) + ":" +
                  boost::lexical_cast<string>(xval[xval.size() - 1]) + "> ";
    }
    // fig_spec += " gnuplot<set key horizontal center>";
    if (expr.size() > 5)
        fig_spec += " gnuplot<set key outside>";
      // fig_spec += " gnuplot<set label 'MS-COOP' at graph 0.5,0.5 center font
      // 'Arial,26'\n set ytics 0.2\n set key horizontal above>";
      //fig_spec +=" gnuplot<set label 'MS-COOP' at graph 0.5,0.5 center font 'Arial,26'\n set ytics 0.2\n unset key>";


    fig_spec += " lw<7>";
    fig_spec += " color<YiZhang16,16>";
    string lege;
    for (size_t i = 0; i < expr.size(); i++) lege += expr[i] + ";";
    fig_spec += " legend<" + lege + ">";
    plot(fig_spec, line_series(xval, lines));
    return;
  }
}

sp_cx_mat coop_grape::update_rf_ham(Eigen::Map<const mat> &m, int index,
                                    size_t step, size_t channel,
                                    size_t nchannels, double kx,
                                    double ky) const {
  double ux = m(2 * nchannels * step + 2 * channel, index);
  double uy = m(2 * nchannels * step + 2 * channel + 1, index);
  // Rf inhom
  ux *= kx;
  uy *= ky;
  return ux * superop_.rf_ctrl.Lx[channel] + uy * superop_.rf_ctrl.Ly[channel];
}

sp_cx_mat coop_grape::update_rf_ham(Eigen::Map<const vec> &v, size_t step,
                                    size_t channel, size_t nchannels, double kx,
                                    double ky) const {
  //double ux = v(2 * nchannels * step + 2 * channel);
  //double uy = v(2 * nchannels * step + 2 * channel + 1);
  //// Rf inhom
  //ux *= kx;
  //uy *= ky;
  //return ux * superop_.rf_ctrl.Lx[channel] + uy * superop_.rf_ctrl.Ly[channel];
  double ux = v(step);
  // Rf inhom
  ux *= kx;
  return ux * superop_.rf_ctrl.Lx[channel];
}

vec spec_avg(const sol::table &t) {
  vec v_avg;
  for (size_t i = 0; i < t.size(); i++) {
    sol::object val = t[i + 1];
    vec v = val.as<vec>();
    if (i == 0) {
      v_avg = v;
      continue;
    }
    v_avg += v;
  }

  v_avg /= (double)t.size();

  return v_avg;
}
vec abs(const vec &a, const vec &b) {
  cx_vec c(a.size());
  c.real() = a;
  c.imag() = b;
  return c.cwiseAbs();
}
}  // namespace oc
}  // namespace ssl
