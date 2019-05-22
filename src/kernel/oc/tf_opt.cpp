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

// Created by  Cecilia Curreli, Technical University Munich on 18-7-4.
#ifdef TENSORFLOW_ENABLED
#include <Python.h>
#include "tensorflow/core/public/session.h"
#include "tensorflow/cc/ops/standard_ops.h" //maybe we don't need this

#include <matrix_exp_op/kernels/matrix_exp_kernel.cc>
#include <matrix_exp_op/kernels/matrix_exp_op.cc>
#include <matrix_exp_op/kernels/matrix_exp_op.h>
#include <matrix_exp_op/ops/matrix_exp.cc>
//#include <tuple>

#include "tf_opt.h"
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

typedef std::chrono::high_resolution_clock::time_point TimeVar;
#define duration(a) std::chrono::duration_cast<std::chrono::seconds>(a).count()
#define timeNow() std::chrono::high_resolution_clock::now()

namespace ssl {
namespace oc {

tf_opt::tf_opt(spin_system &sys)
    : rf_(nullptr) {

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

tf_opt::~tf_opt() {
}

void tf_opt::create_graph(std::string g_type) {

  int n_channels = sys_->nspins();
  const int size = (int) (init_state_).size();
  std::string graph_definition = "opt_graph.pb";
  bool measure_time = true;
  TimeVar start;

  if (measure_time) { start = timeNow(); }

  bool use_gpu = 0;

  //call system which will call python function to construct the graph
  if (call_system(size, n_channels, dt, graph_definition, g_type)) {
    std::cout << "\n\n" << std::endl;
    string s1 = str(boost::format("Time needed for graph computation: [%.3f] seconds.\n\n") %
        ((duration(timeNow() - start))));
    ssl_color_text("info", s1);
  }
  if (measure_time) { start = timeNow(); }

  //load the graph
  TF_CHECK_OK(ReadBinaryProto(Env::Default(), graph_definition, &graph_def));

  //assign Device after defining the graph.
  if (opts.target.empty() && use_gpu) {
    std::cout << "empty" << std::endl;
    graph::SetDefaultDevice(use_gpu ? "/gpu:0" : "/cpu:0", &graph_def);
  }

  // create a new session
  TF_CHECK_OK(NewSession(opts, &session));
  // Load graph into session
  TF_CHECK_OK(session->Create(graph_def));

  // this can be used to verify which devices are found
  /*
  std::vector<DeviceAttributes> resp;
  TF_CHECK_OK(session->ListDevices(&resp));
  bool has_gpu = false;
  for (const auto& dev : resp) {
      if (dev.device_type() == "GPU") { has_gpu = true; std::cout << "Device found: GPU " << std::endl;}
      if (dev.device_type() == "CPU") { std::cout << "Device found: CPU " << std::endl; }
  }
  */
  //

  create_tf_variables_list();
  initialize_placeholder();

  if (measure_time) {
    string
        s1 = str(boost::format("Time needed to load the graph: [%.3f] seconds.\n\n") % ((duration(timeNow() - start))));
    ssl_color_text("info", s1);
  }

}

void tf_opt::update_tf_variables(const double *x) {

  int n_channels = sys_->nspins();
  Tensor *u_pointer = nullptr;

  for (int i = 0; i < nsteps; ++i) {
    for (int c = 0; c < n_channels; ++c) {

      u_pointer = &std::get<1>(variables_list[i * n_channels * 2 + c * 2]);
      u_pointer->scalar<double>()() = x[i * n_channels * 2 + c * 2];
      //memcpy could be a bit faster.
      //memcpy(u_pointer->flat<double>().data(),
      // &raw_val[i * n_channels * 2 + c * 2], sizeof(double)); //u_pointer->scalar<double>()() = -250*2*M_PI ;

      u_pointer = &std::get<1>(variables_list[i * n_channels * 2 + c * 2 + 1]);
      u_pointer->scalar<double>()() = x[i * n_channels * 2 + c * 2 + 1];
      //memcpy could be a bit faster.
      //memcpy(u_pointer->flat<double>().data(),
      // &raw_val[i * n_channels * 2 + c * 2 + 1], sizeof(double)); // u_pointer->scalar<double>()() = 0.0 ;

    }
  }

  TF_CHECK_OK(session->Run(variables_list, {}, {"init_all_vars_op"}, nullptr));

}

double tf_opt::objfunc(const vector<double> &x, vector<double> &grad) {

  rf_->update_raw_data(x.data());
  //update variables in graph with values of x
  update_tf_variables(x.data());

  //get gradients
  TF_CHECK_OK(session->Run(placeholder_list, gradient_names, {}, &outputs)); //add "trace_real" to the gradient_names
  //save gradient in grad
  Tensor *u_pointer = nullptr;
  for (int i = 0; i < (outputs.size() - 1); ++i) {
    u_pointer = &outputs[i];
    grad[i] = -u_pointer->scalar<double>()();
  }
  double fidelity = (outputs[outputs.size() - 1].scalar<double>()());// - 1;
  outputs.clear();
  //TO DO
//            TF_CHECK_OK(session->Run(placeholder_list, {"fidelity", "psi" + std::to_string(nsteps), "Ux0_0","Uy0_0", "trace_real", "trace"}, {},  &outputs));
////            double fidelity = (outputs[0].scalar<double>()());
//            std::cout << "Ux in x:    " << x[0] << ";    Ux in tensorflow:    " << outputs[2].scalar<double>()()<<std::endl;
//            std::cout << "Uy in x:     " << x[0 + 1] << ";    Uy in tensorflow:    " << outputs[3].scalar<double>()()<< std::endl;
//                std::cout << "grad wrt to first Ux:    " << grad[0] << std::endl;
//                std::cout << "grad wrt to first Uy:     " << grad[0 + 1] << std::endl;
//    //            std::cout << "trace.real():    " << outputs[0].scalar<double>()()<< std::endl;
//            std::cout << "fidelity tf:    " << (outputs[0].scalar<double>()()) << std::endl;
//            std::cout << "fidelity:    " << fidelity << std::endl;
//            std::cout << "trace_real tf:    " << (outputs[4].flat<double>()) << std::endl;
//            std::cout << "trace tf:    " << (outputs[5].flat<complex<double>>()) << std::endl;
//            outputs.clear();


//            std::cout << "final state :   " <<outputs[1].flat<complex<double>>() << std::endl;
//            std::cout << "target state :   " <<Eigen::VectorXcd(targ_state_) << std::endl;


  std::cout << boost::format("==> %04d  [%.4f]\n") % (++opt_.iter_count) % fidelity;
  opt_.vf.push_back(fidelity);
  return fidelity;

}


//        void tf_opt::evaluation(const sol::table &t) {
//
//            double width = retrieve_table_double("width", t); // unit in ms.
//            nsteps = (size_t) (retrieve_table_double("step", t));
//            dt = width * 1e-3 / double(nsteps); // into s.
//            if  (Eigen::MatrixXcd(superop_.R).isZero()) {relax = false;}
//            else {relax = true;}
//
//
//
//            bool measure_time = 1;
//            TimeVar start;
//            if(measure_time) { start=timeNow(); }
//
//            assign_state(t);
//            assign_nlopt(t);
//            assign_pulse(t);
//            assign_aux_var();
//
//            create_graph("evaluation");
//
//            rf_->convert2(_ux_uy);
//            vector<double> x = rf_->clone_raw_data();
//
//            //create vector with the names of the grafients to fetch gradient with a call to tensroflow run
//            int n_channels = sys_->nspins();
//            gradient_names = std::vector<string> (nsteps*n_channels*2 + 1); //here we have all nodes we want to fetch for the optimization
//            string s;
//            for (int i = 0; i < nsteps; ++i) {
//                for (int c = 0; c < n_channels; ++c) {
//                    s = std::to_string(c) + "_" + std::to_string(i);
//                    gradient_names[i*n_channels*2 + c*2] = "Cecilia_gradient/Complex_Ux" + s + "_grad/Reshape";
//                    gradient_names[i*n_channels*2 + c*2 + 1] = "Cecilia_gradient/Complex_Uy" + s + "_grad/Reshape";
//                }
//            }
//            gradient_names[nsteps*n_channels*2] = "trace_real";
//
//            // optimization processing.
//            nlopt::opt opt(opt_.algo, rf_->get_dims());
//            opt.set_xtol_rel(opt_.tol_f);
//            if (opt_.max_eval > 0)
//                opt.set_maxeval(opt_.max_eval);
//            if (opt_.max_time > 0)
//                opt.set_maxtime(opt_.max_time);
//            if (opt_.stopval > 0)
//                opt.set_stopval(opt_.stopval);
//
//            opt.set_max_objective(objfunc, this);
//
//            // boundary limitation.
//            // opt_amplitude_constraint(opt);
//
//            double max_f;
//            nlopt::result result = opt.optimize(x, max_f);
//
//            if (result == nlopt::SUCCESS)
//                ssl_color_text("info", "pulse optimization succeed.\n");
//            if(measure_time) {
//                string s1 = str(boost::format("Time needed for the Tensorflow optimization process is: [%.3f] seconds.\n") % duration(timeNow()-start));
//                ssl_color_text("info", s1);
//            }
//            if (result == nlopt::FAILURE)
//                ssl_color_text("err", "pulse optimization failed.\n");
//            if (result == nlopt::MAXEVAL_REACHED)
//                ssl_color_text("warn", "pulse optimization terminated due to maximum iterations limitation.\n");
//
//            opt_.max_f = max_f;
//            rf_->update_raw_data(x.data());
//            update_tf_variables(x.data());
//
//            tf_opt::h5write();
//
//            //session->Close(); //this eliminates all the variables
//            //delete session;
//
//            return *rf_;
//        }

seq_block &tf_opt::optimize(const sol::table &t) {

  double width = retrieve_table_double("width", t); // unit in ms.
  nsteps = (size_t) (retrieve_table_double("step", t));
  dt = width * 1e-3 / double(nsteps); // into s.
  if (Eigen::MatrixXcd(superop_.R).isZero()) { relax = false; }
  else { relax = true; }

  bool measure_time = 1;
  TimeVar start;
  if (measure_time) { start = timeNow(); }

  assign_state(t);
  assign_nlopt(t);
  assign_pulse(t);
  assign_aux_var();

  create_graph("optimize");

  rf_->convert2(_ux_uy);
  vector<double> x = rf_->clone_raw_data();

  //create vector with the names of the grafients to fetch gradient with a call to tensroflow run
  int n_channels = sys_->nspins();
  gradient_names =
      std::vector<string>(nsteps * n_channels * 2 + 1); //here we have all nodes we want to fetch for the optimization
  string s;
  for (int i = 0; i < nsteps; ++i) {
    for (int c = 0; c < n_channels; ++c) {
      s = std::to_string(c) + "_" + std::to_string(i);
      gradient_names[i * n_channels * 2 + c * 2] = "Cecilia_gradient/Complex_Ux" + s + "_grad/Reshape";
      gradient_names[i * n_channels * 2 + c * 2 + 1] = "Cecilia_gradient/Complex_Uy" + s + "_grad/Reshape";
    }
  }
  gradient_names[nsteps * n_channels * 2] = "trace_real";

  // optimization processing.
  nlopt::opt opt(opt_.algo, rf_->get_dims());
  opt.set_xtol_rel(opt_.tol_f);
  if (opt_.max_eval > 0)
    opt.set_maxeval(opt_.max_eval);
  if (opt_.max_time > 0)
    opt.set_maxtime(opt_.max_time);
  if (opt_.stopval > 0)
    opt.set_stopval(opt_.stopval);

  opt.set_min_objective(objfunc, this);

  // boundary limitation.
  // opt_amplitude_constraint(opt);

  double max_f;
  nlopt::result result = opt.optimize(x, max_f);

  if (result == nlopt::SUCCESS)
    ssl_color_text("info", "pulse optimization succeed.\n");
  if (measure_time) {
    string s1 = str(boost::format("Time needed for the Tensorflow optimization process is: [%.3f] seconds.\n")
                        % duration(timeNow() - start));
    ssl_color_text("info", s1);
  }
  if (result == nlopt::FAILURE)
    ssl_color_text("err", "pulse optimization failed.\n");
  if (result == nlopt::MAXEVAL_REACHED)
    ssl_color_text("warn", "pulse optimization terminated due to maximum iterations limitation.\n");

  opt_.max_f = max_f;
  rf_->update_raw_data(x.data());
  update_tf_variables(x.data());

  tf_opt::h5write();

  //session->Close(); //this eliminates all the variables
  //delete session;

  return *rf_;
}

void tf_opt::projection(const sol::table &t) { //results are not correct
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
}

void tf_opt::create_tf_variables_list() {

  int n_channels = sys_->nspins();
  variables_list = std::vector<std::pair<std::string, Tensor> >(nsteps * n_channels * 2);

  Tensor *u_pointer = nullptr;
  for (int i = 0; i < nsteps; ++i) {
    for (int c = 0; c < n_channels; ++c) {

      u_pointer = new Tensor(DT_DOUBLE, tensorflow::TensorShape({1}));
      variables_list[i * n_channels * 2 + c * 2] =
          {"Assignx" + std::to_string(c) + "_" + std::to_string(i), *u_pointer};

      u_pointer = new Tensor(DT_DOUBLE, tensorflow::TensorShape({1}));
      variables_list[i * n_channels * 2 + c * 2 + 1] =
          {"Assigny" + std::to_string(c) + "_" + std::to_string(i), *u_pointer};

    }
  }

}

void tf_opt::initialize_placeholder() {

  const int size = (int) init_state_.size();
  int n_channels = sys_->nspins();

  placeholder_list = std::vector<std::pair<std::string, Tensor> >(2 * n_channels + 2 + 1);

  Tensor *Lx = nullptr;
  Tensor *Ly = nullptr;

  //get values for placeholder by copying into tensors
  Tensor psi0(DT_COMPLEX128, tensorflow::TensorShape({size, 1}));
  auto ps = Eigen::VectorXcd(init_state_);
  memcpy(psi0.flat<complex<double>>().data(), ps.data(), sizeof(complex<double>) * size);
  //std::cout << "initial state:\n" << ps << std::endl;

  Tensor target(DT_COMPLEX128, tensorflow::TensorShape({size, 1}));
  auto psn = Eigen::VectorXcd(targ_state_);
  memcpy(target.flat<complex<double>>().data(), psn.data(), sizeof(complex<double>) * size);
  //std::cout << "target state:\n" << psn << std::endl;

  Tensor L0(DT_COMPLEX128, tensorflow::TensorShape({size, size}));
  auto L_free = Eigen::MatrixXcd(superop_.L0);
  //this works only if L0 is always symmetric! otherwise we would copy a column mayor matrix into a row major matrix
  //memcpy((L0.flat<complex<double>>()).data(), L_free.data(), sizeof(complex<double>) * size * size);
  auto mapL0 = Eigen::Map<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor  >>(
      (L0.flat<complex<double>>()).data(), size, size); /* tensorflow::Tensor is always row-major */
  mapL0 = L_free;

  placeholder_list[n_channels * 2] = {"psi0", psi0};
  placeholder_list[n_channels * 2 + 1] = {"target", target};
  placeholder_list[n_channels * 2 + 2] = {"L0", L0};

  if (relax) {
    Tensor R(DT_COMPLEX128, tensorflow::TensorShape({size, size}));
    auto R_sys = Eigen::MatrixXcd(superop_.R);

    auto mapR = Eigen::Map<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor  >>(
        (R.flat<complex<double>>()).data(), size, size); /* tensorflow::Tensor is always row-major */
    mapR = R_sys;
    placeholder_list.push_back({"R", R});
    std::cout << "Ro matrix:\n" << superop_.R << std::endl;
    std::cout << "R0 matrix:\n" << Eigen::MatrixXcd(superop_.R) << std::endl;
    std::cout << "R0 tensor:\n" << R.matrix<complex<double>>() << std::endl;

  }

  for (int c = 0; c < n_channels; ++c) {
    Lx = new Tensor(DT_COMPLEX128, tensorflow::TensorShape({size, size}));
    Ly = new Tensor(DT_COMPLEX128, tensorflow::TensorShape({size, size}));

    auto mapLx = Eigen::Map<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor  >>(
        (Lx->flat<complex<double>>()).data(), size, size); /* tensorflow::Tensor is always row-major */
    auto mapLy = Eigen::Map<Eigen::Matrix<complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor  >>(
        (Ly->flat<complex<double>>()).data(), size, size); /* tensorflow::Tensor is always row-major */

    //solution 1, maybe not so fast
    mapLx = Eigen::MatrixXcd(superop_.rf_ctrl.Lx[c]);
    mapLy = Eigen::MatrixXcd(superop_.rf_ctrl.Ly[c]);

    //solution 2, could be faster.
    //memcpy(Ly->matrix<complex<double>>().data(),
    //       (Eigen::MatrixXcd(superop_.rf_ctrl.Ly[c])).data(), sizeof(complex<double>) * size * size);
    //mapLy.transposeInPlace();
    //memcpy(Lx->matrix<complex<double>>().data(),
    //       (Eigen::MatrixXcd(superop_.rf_ctrl.Lx[c])).data(), sizeof(complex<double>) * size * size);
    //mapLx.transposeInPlace();

    placeholder_list[c] = {"Lx" + std::to_string(c), *Lx};
    placeholder_list[c + n_channels] = {"Ly" + std::to_string(c), *Ly};

  }
}

bool tf_opt::call_system(int size,
                         int n_channels,
                         double dt,
                         std::string graph_definition,
                         std::string graph_type) //const int* inputs
{
  if (system(nullptr)) {
    string s = "cd " + g_project_path + "; cp -r share/spin-scenario/tf_files" + " " + g_spin_scenario + "; cd "
        + g_spin_scenario;
    int res = system(s.c_str());
    string command = "cd tf_files; python3 -c 'import graph_generator as gg; print (gg.compute_graph_";
    command += graph_type + "(" + to_string(size) + ", " + to_string(nsteps) + ", " + to_string(n_channels) + ", " +
        to_string(dt) + ", " + to_string(relax) + "))'";
    res = system(command.c_str());
    struct stat buf;
    std::cout << "output of system " << res << std::endl; //= 256 if error //= 0 if everything ok
    if (!(stat(graph_definition.c_str(), &buf) == 0)) {
      string s = "\n           Error in system() call to python code. Graph was not created!!";
      throw std::runtime_error(s.c_str());
      return 0;
    } else {
      return 1;
    }
  } else {
    string s = "!!command processor doesn't exists!! Failed call to system()!!";
    throw std::runtime_error(s.c_str());
    return 0;
  }
}
double tf_opt::maxf() const {
  return opt_.max_f;
}

double tf_opt::objfunc(const vector<double> &x, vector<double> &grad, void *func) {
  return ((tf_opt *) func)->objfunc(x, grad);
}

//double grape::objfunc_broadband(const vector<double> &x, vector<double> &grad, void *func) {
//    return ((grape *) func)->objfunc_broadband(x, grad);
//}

void tf_opt::assign_state(const sol::table &t) {
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

/*  if (obj.get_type() == sol::type::table) {
    sol::table targ_table = obj.as<sol::table>();
    for (size_t i = 0; i < targ_table.size(); i++) {
      sol::object val = targ_table[i + 1];
      sp_cx_vec targ = val.as<sp_cx_vec>();
      targ = levante_ernst(targ);
      targ_list_.push_back(targ);
    }
  }*/
}

void tf_opt::assign_nlopt(const sol::table &t) {
  opt_.algo = nlopt::LD_LBFGS; // default algorithm.
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

void tf_opt::assign_pulse(const sol::table &t) {
  double width = retrieve_table_double("width", t); // unit in ms.
  size_t nsteps = (size_t) (retrieve_table_double("step", t));
  string pattern = "random_spline";
  if (is_retrievable("init_pattern", t))
    pattern = retrieve_table_str("init_pattern", t);

  double dt = width * 1e-3 / double(nsteps); // into s.

  string str_chs = boost::algorithm::join(superop_.rf_ctrl.chs, " ");
  string code = "user_rf = shapedRF{name = 'juncy', width = " + boost::lexical_cast<string>(width) +
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
//    double vals [100];
//    for (int i = 0; i < 100; ++i) {
//        if (i%2==0){vals[i] = -250*2*M_PI;}
//        else{ vals[i] =0.0;}
//    }
//    rf_->update_raw_data(vals);

}

void tf_opt::h5write(string file_name) const {
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

void tf_opt::assign_aux_var() {
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

//void tf_opt::opt_amplitude_constraint(nlopt::opt &opt) {
//    if (opt_.max_amp <= 0)
//        return;
//
//    size_t dim = rf_->get_dims();
//    size_t nsteps = rf_->get_steps();
//    size_t nchannels = rf_->get_channels();
//    vector<double> up_bound(dim);
//    vector<double> low_bound(dim);
//
//    switch (opt_.algo) {
//        case nlopt::LD_MMA:
//        case nlopt::LD_LBFGS: {
//            //if(rf_mode_ ==_amp_phase)
//            //      for (size_t i = 0; i < opt_rf_.nsteps; i++) {
//            //          for (size_t j = 0; j < opt_rf_.nchannels; j++) {
//            //              // amp.
//            //              up_bound[2 * opt_rf_.nchannels * i + 2 * j] = opt_rf_.max_amp * 2 * _pi;
//            //              low_bound[2 * opt_rf_.nchannels * i + 2 * j] = 0;
//
//            //              // phase.
//            //              up_bound[2 * opt_rf_.nchannels * i + 2 * j + 1] = 2 * _pi;
//            //              low_bound[2 * opt_rf_.nchannels * i + 2 * j + 1] = 0;
//            //          }
//            //      }
//
//            if (rf_->mode() == _ux_uy)
//                for (size_t i = 0; i < nsteps; i++) {
//                    for (size_t j = 0; j < nchannels; j++) {
//                        // ux.
//                        up_bound[2 * nchannels * i + 2 * j] = opt_.max_amp * 2 * _pi;
//                        low_bound[2 * nchannels * i + 2 * j] = -opt_.max_amp * 2 * _pi;
//
//                        // uy.
//                        up_bound[2 * nchannels * i + 2 * j + 1] = opt_.max_amp * 2 * _pi;
//                        low_bound[2 * nchannels * i + 2 * j + 1] = -opt_.max_amp * 2 * _pi;
//                    }
//                }
//
//            opt.set_upper_bounds(up_bound);
//            opt.set_lower_bounds(low_bound);
//        }
//            break;
//            /*case nlopt::LD_MMA: {
//                amp_constraint_data data = { 0, 0, opt_rf_.nchannels, opt_rf_.max_amp * 2 * _pi };
//                for (size_t i = 0; i < opt_rf_.nsteps; i++) {
//                    for (size_t j = 0; j < opt_rf_.nchannels; j++) {
//                        data.i = i;
//                        data.j = j;
//                        opt.add_inequality_constraint(amplitude_constraint, &data, opt_.opt_f);
//                    }
//                }
//            }
//            break;*/
//        default:break;
//    }
//}


void tf_opt::cartesian2polar(double amp, double phase, double gx, double gy, double &g_amp, double &g_phase) {
  double ux = amp * cos(phase);
  double uy = amp * sin(phase);
  g_amp = gx * ux + gy * uy;
  g_amp /= amp;
  g_phase = -gx * uy + gy * ux;
}

//double tf_opt::amplitude_constraint(unsigned n, const double *x, double *grad, void *data) {
//    amp_constraint_data *d = reinterpret_cast<amp_constraint_data *>(data);
//
//    double ux = x[2 * d->chs * d->i + 2 * d->j];
//    double uy = x[2 * d->chs * d->i + 2 * d->j + 1];
//
//    return ux * ux + uy * uy - d->amp * d->amp;
//}


//double tf_opt::objfunc_broadband(const vector<double> &x, vector<double> &grad) {
//    rf_->update_raw_data(x.data());
//    int N = superop_.L0s.size();
//    vec phi = vec::Zero(N);
//
//    size_t nsteps = rf_->get_steps();
//    size_t nchannels = rf_->get_channels();
//    double dt = rf_->width_in_ms() * 1e-3 / double(nsteps); // into s.
//
//#pragma omp parallel for
//    for (int p = 0; p < N; p++) {
//        int id = omp_get_thread_num();
//        traj_omp_[id].forward[0] = init_state_;
//        traj_omp_[id].backward[nsteps] = targ_list_[p];
//        //cout << id << "\n";
//        sp_cx_mat L;
//        sp_cx_mat L0 = superop_.L0s[p] + ci * superop_.R;
//        double kx = 1, ky = 1;
//        sp_cx_vec rho = traj_omp_[id].forward[0];
//        for (size_t i = 0; i < nsteps; i++) {
//            L = L0;
//            for (size_t j = 0; j < nchannels; j++)
//                L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
//            traj_omp_[id].forward[i + 1] = ssl::spinsys::step(rho, L, dt);
//            rho = traj_omp_[id].forward[i + 1];
//        }
//        rho = traj_omp_[id].backward[nsteps];
//        for (int i = nsteps - 1; i >= 0; i--) {
//            L = L0.adjoint();
//            for (size_t j = 0; j < nchannels; j++)
//                L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
//            traj_omp_[id].backward[i] = ssl::spinsys::step(rho, L, -dt);
//            rho = traj_omp_[id].backward[i];
//        }
//        sp_cx_mat Gx, Gy, tmp;
//        int k = 0;
//        for (size_t i = 0; i < nsteps; i++) {
//            L = L0;
//            for (size_t j = 0; j < nchannels; j++)
//                L += update_rf_ham(x.data(), i, j, nchannels, kx, ky);
//            tmp = traj_omp_[id].backward[i + 1].adjoint() * ssl::spinsys::propagator(L, dt);
//            for (size_t j = 0; j < nchannels; j++) {
//                Gx = propagator_derivative(L, superop_.rf_ctrl.Lx[j], dt);
//                Gy = propagator_derivative(L, superop_.rf_ctrl.Ly[j], dt);
//
//                superop_.grad_bb[p][k] = traced(tmp * Gx * traj_omp_[id].forward[i]).real();
//                superop_.grad_bb[p][k + 1] = traced(tmp * Gy * traj_omp_[id].forward[i]).real();
//                k += 2;
//            }
//        }
//        phi[p] = transfer_fidelity(traj_omp_[id].forward[nsteps], targ_list_[p]);
//        /// transfer_fidelity(targ_list_[p], targ_list_[p]);
//    } // end parallel for.
//
//    double val = phi.sum() / (double) N;
//    double alpha = 0;
//
//    int k = 0;
//    for (size_t i = 0; i < nsteps; i++) {
//        for (size_t j = 0; j < nchannels; j++) {
//            double gx = 0;
//            double gy = 0;
//            for (int p = 0; p < N; p++) {
//                gx += superop_.grad_bb[p][k];
//                gy += superop_.grad_bb[p][k + 1];
//            }
//            grad[k] = gx / (double) N;
//            grad[k + 1] = gy / (double) N;
//
//            //double gxx = grad[k];
//            //double gyy = grad[k + 1];
//
//            //if (opt_rf_.rf_mode == _amp_phase)  // Translate the control gradient into phase gradient
//            //cartesian2polar(x[k], x[k + 1], gxx, gyy, grad[k], grad[k + 1]);
//
//            // rf power reduction.
//            /*double ux, uy;
//            ux = x[2 * opt_rf_.nchannels * i + 2 * j];
//            uy = x[2 * opt_rf_.nchannels * i + 2 * j + 1];
//            grad[k] += 2.0 * alpha*ux*opt_rf_.dt;
//            grad[k+1] += 2.0 * alpha*uy*opt_rf_.dt;*/
//
//            k += 2;
//        }
//    }
//
//    double phi_rf = rf_->average_power();
//    double PHI1 = val;
//    double PHI2 = alpha * phi_rf * dt;
//    std::cout << boost::format("==> %04d  [%.4f] [%.4f]\n") % (++opt_.iter_count) % PHI1 % PHI2;
//    opt_.vf.push_back(PHI1 - PHI2);
//    return (PHI1 - PHI2);
//}
sp_cx_mat tf_opt::update_rf_ham(const double *x,
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
}

}
#endif  // TENSORFLOW_ENABLED
