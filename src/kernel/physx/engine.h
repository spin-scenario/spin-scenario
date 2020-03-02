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

#include <kernel/spinsys/spin_system.h>
using namespace ssl::spinsys;
#include <kernel/sample/phantom.h>
using namespace ssl::sample;

//#define DENSE_MATRIX_COMPUTE 1
#define ARRAYFIRE_USE 1

#ifdef ARRAYFIRE_USE
#include <af/util.h>
#include <arrayfire.h>
#endif

using namespace ssl::utility;

namespace ssl {
namespace physx {

enum physx_model {
  _quantum_cpu,
  _quantum_gpu,
  _bloch_gpu
};

struct unified_spinsys {
  const spin_system *p_sys;
  sp_cx_mat L0; // free ham.
  sp_cx_mat Lz0;
  sp_cx_mat R; // relaxation op.
  rf_ham rf_ctrl; // rf hams.
  sp_cx_vec rho0;
  double B0; // unit in Tesla.
#ifdef DENSE_MATRIX_COMPUTE
  cx_vec det;
#else
  sp_cx_vec det;
#endif

  void init(const spin_system &sys) {
    p_sys = &sys;
    B0 = sys.get_magnet_field();
    rf_ctrl = sys.rf_hamiltonian();
    for (size_t i = 0; i < rf_ctrl.channels; i++) {
      rf_ctrl.Lx_dense[i] = rf_ctrl.Lx[i].toDense();
      rf_ctrl.Ly_dense[i] = rf_ctrl.Ly[i].toDense();
    }

    L0 = sys.free_hamiltonian();
    R = ci * sys.relaxation();

    string acq_ch = g_seq_param->acq_channel;
    //rho0 = sys.state(acq_ch + " Iz");
    rho0 = sys.equilibrium_state();
    Lz0 = sys.op(acq_ch + " Iz").cast<cd>();

#ifdef DENSE_MATRIX_COMPUTE
    det = sys.state("1H I+").toDense(); // temporarily used.
    //det = sys.smart_state("I1y").toDense(); // temporarily used.

#else
if(g_seq_param->observer.empty())
    det = sys.state(acq_ch + " I+"); // temporarily used.
else
{
    cout<<"###\n";
    det= sys.smart_state(g_seq_param->observer);

}
#endif
    double norm_rho0 = rho0.norm();
    double norm_det =det.norm();
    det *= norm_rho0 / norm_det; // this make norm(det)==norm(rho0), please DO NOT change rho0.
    levante_ernst_correction(rho0);
    levante_ernst_correction(det);
  }
};

struct each_spinsys {
#ifdef DENSE_MATRIX_COMPUTE
  cx_mat L;
  cx_mat L0; // free ham of this isochromat.
  cx_mat Lz0;
  cx_mat R; // relaxation op.
  cx_vec rho;
#else
  sp_cx_mat L;
  sp_cx_mat L0; // free ham of this isochromat.
  sp_cx_mat Lz0;
  sp_cx_mat R; // relaxation op.
  sp_cx_vec rho;
#endif
  vec3 pos;
  double pd;
  double dB;
  cx_vec sig;
};

struct rf_const {
  vec2 u = vec2::Zero(); // x/y components.
  double df = 0; // Hz
  string channel = "";
};

struct grad_const {
  vec3 v = vec3::Zero(); // x/y/z components.
};

struct acq_const {
  bool adc = false;
  bool last = false;
  int nps = 0;
  int index = -1;
  bool shot = false;
};

struct seq_const {
  vector<rf_const> rf;
  grad_const grad;
  acq_const acq;
  bool delay_if = false;
  bool rf_if = false;
};

struct seq_step {
  seq_const cur;
  double dt;
};

class engine {
 public:
  engine(const sol::table &t);
  ~engine();
  void evolution(timeline dt, const seq_const &ctrl);
  //void set_observer(const sp_cx_vec &rho);
  //void set_observer(const string expr);
  sol::object process_signal();

 private:
  void init_ensemble(const phantom *p_phantom);
  void TFInitEnsemble(const phantom *p_phantom);
  void evolution_for_each(double dt, const seq_const &ctrl, each_spinsys &each);
  cx_vec accu_signal();
  cx_vec accu_signal_bloch_gpu();
  void gpu_signal_transfer(double *fidx, double *fidy);
 private:
  physx_model physx_model_;
  phantom *p_phantom_;
  unified_spinsys unified_spinsys_;
  vector<each_spinsys> ensemble_;

  vector<seq_step> seq_step_list_;

  vector<cx_vec> raw_signal_;
  vector<cx_vec> raw_signal_tmp_;

  int nRxs_;
  int nTxs_;

  #ifdef ARRAYFIRE_USE
  af::array gPos_;
  af::array gdB0s_;

  af::array tx_sens_amp_;
  af::array tx_sens_phi_;
  af::array rx_sens_amp_;
  af::array rx_sens_phi_;

  af::array rx_sens_cx_;
  af::array rx_sens_cy_;

  af::array r1_;
  af::array r2_;
  af::array pd_;
  af::array ux_;
  af::array uy_;

  vector<double> raw_r1_;
  vector<double> raw_r2_;
  vector<double> raw_rho_;
  vector<double> raw_loc_;
  vector<double> raw_dB0_;
  vector<double> raw_tx_B1_amp_;
  vector<double> raw_tx_B1_phase_;
  vector<double> raw_rx_B1_amp_;
  vector<double> raw_rx_B1_phase_;



  af::array mx_;
  af::array my_;
  af::array mz_;

  //cx_mat signal_;
  //vector<cx_mat> signals_;
#endif
};

extern engine *g_engine;
}
}
