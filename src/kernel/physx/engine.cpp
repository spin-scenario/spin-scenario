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

#include "engine.h"
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;

namespace ssl {
namespace physx {

engine *g_engine = nullptr;
engine::engine(const sol::table &t) {
  string acc_str = "cpu";  // default option.
  if (is_retrievable("acc", t)) {
    acc_str = retrieve_table_str("acc", t);
  }
  boost::to_lower(acc_str);
  if (acc_str == "gpu")
    physx_model_ = _bloch_gpu;  // currently we support GPU version only for
                                // Bloch model (via Arrayfire lib).
  else if (acc_str == "cpu") {
    physx_model_ =
        _quantum_cpu;  // for the cpu option, we default use quantum model.
  }

  if (is_retrievable("spinsys", t)) {
    // general one spin system case.
    sol::object par = retrieve_table("spinsys", t);
    const spin_system &par_sys = par.as<const spin_system &>();
    unified_spinsys_.init(par_sys);
  } else {
    // g_lua->script("_sys = spin_system{B0 = '3 T', spins = '1H'}");
    g_lua->script("_sys = spin_system{spin = '1H'}");
    unified_spinsys_.init((*g_lua)["_sys"]);
  }

  p_phantom_ = nullptr;
  if (is_retrievable("phantom", t)) {
    string par = retrieve_table_str("phantom", t);
    p_phantom_ = new phantom(par.c_str());
  }
  init_ensemble(p_phantom_);
  delete[] p_phantom_;
}

engine::~engine() {}

void engine::init_ensemble(const phantom *p_phantom) {
  double ppm = 0;  // 0.05;
  if (physx_model_ == _quantum_cpu) {
    if (!p_phantom) {
      ssl_color_text("warn",
                     "no phantom specified, ONLY one spin system used.\n");
      each_spinsys each;
      ensemble_ = vector<each_spinsys>(1, each);
#ifdef DENSE_MATRIX_COMPUTE
      ensemble_[0].pos = vec3::Zero();
      ensemble_[0].rho = unified_spinsys_.rho0.toDense();
      ensemble_[0].L0 =
          unified_spinsys_.L0
              .toDense();  //////////////////////////////////////////////////////////////////////////
      ensemble_[0].Lz0 = unified_spinsys_.Lz0.toDense();
      ensemble_[0].R = unified_spinsys_.R.toDense();
      ensemble_[0].pd = 1;
#else
      ensemble_[0].pos = vec3::Zero();
      ensemble_[0].rho = unified_spinsys_.rho0;
      ensemble_[0].L0 =
          unified_spinsys_
              .L0;  //////////////////////////////////////////////////////////////////////////
      ensemble_[0].Lz0 = unified_spinsys_.Lz0;
      ensemble_[0].R = unified_spinsys_.R;
      ensemble_[0].pd = 1;
#endif
      return;
    }

    int num = p_phantom->isochromats_.size();
    each_spinsys each;
    ensemble_ = vector<each_spinsys>(num, each);
    omp_set_num_threads(omp_core_num);
#pragma omp parallel for
    for (int i = 0; i < num; i++) {
      // int id = omp_get_thread_num();
#ifdef DENSE_MATRIX_COMPUTE
      ensemble_[i].pos = p_phantom->isochromats_[i].position();
      ensemble_[i].rho = unified_spinsys_.rho0.toDense();
      ensemble_[i].L0 = unified_spinsys_.L0.toDense();
      ensemble_[i].Lz0 = unified_spinsys_.Lz0.toDense();
      ensemble_[i].R = unified_spinsys_.R.toDense();
      ensemble_[i].R.setZero();
      ensemble_[i].R(2, 0) = ci * p_phantom->isochromats_[i].data[r1];
      ensemble_[i].R(1, 1) = -ci * p_phantom->isochromats_[i].data[r2];
      ensemble_[i].R(2, 2) = -ci * p_phantom->isochromats_[i].data[r1];
      ensemble_[i].R(3, 3) = -ci * p_phantom->isochromats_[i].data[r2];
      //      cout<<ensemble_[i].R<<"\n";
      //      int a;
      //      cin>>a;
      ensemble_[i].dB = p_phantom->isochromats_[i].data[dB0] * ppm *
                        unified_spinsys_.B0 / 1e3;  // into mT
      ensemble_[i].pd = p_phantom->isochromats_[i].data[pd];
#else
      ensemble_[i].pos = p_phantom->isochromats_[i].position();
      ensemble_[i].rho = unified_spinsys_.rho0;
      ensemble_[i].L0 = unified_spinsys_.L0;
      ensemble_[i].Lz0 = unified_spinsys_.Lz0;
      ensemble_[i].R = unified_spinsys_.R;
      ensemble_[i].R.setZero();
      ensemble_[i].R.coeffRef(2, 0) = ci * p_phantom->isochromats_[i].data[r1];
      ensemble_[i].R.coeffRef(1, 1) = -ci * p_phantom->isochromats_[i].data[r2];
      ensemble_[i].R.coeffRef(2, 2) = -ci * p_phantom->isochromats_[i].data[r1];
      ensemble_[i].R.coeffRef(3, 3) = -ci * p_phantom->isochromats_[i].data[r2];
      ensemble_[i].dB = p_phantom->isochromats_[i].data[dB0] * ppm *
                        unified_spinsys_.B0 / 1e3;  // into mT.
      ensemble_[i].pd = p_phantom->isochromats_[i].data[pd];
#endif
    }
  }
  if (physx_model_ == _bloch_gpu) {
#ifdef ARRAYFIRE_USE
    af::setDevice(0);
    af::info();
    int num = p_phantom->isochromats_.size();
    for (int i = 0; i < num; i++) {
      const isochromat &iso = p_phantom->isochromats_[i];
      raw_r1_.push_back(iso.data[r1]);
      raw_r2_.push_back(iso.data[r2]);
      raw_rho_.push_back(iso.data[pd]);

      vector<double> pos(3, 0);
      pos[0] = iso.data[cx];
      pos[1] = iso.data[cy];
      pos[2] = iso.data[cz];
      std::copy(pos.begin(), pos.end(), std::back_inserter(raw_loc_));

      raw_dB0_.push_back(iso.data[dB0]);

      nTxs_ = 1;  // only for single tx coil.
      nRxs_ = 1;  // only for single rx coil.

      for (size_t k = 0; k < nTxs_; k++) {
        // cd val = iso->sens_tx[k];
        raw_tx_B1_amp_.push_back(1);
        raw_tx_B1_phase_.push_back(0);
      }

      for (size_t k = 0; k < nRxs_; k++) {
        // cd val = iso->sens_rx[k];
        raw_rx_B1_amp_.push_back(1);
        raw_rx_B1_phase_.push_back(0);
      }
    }
    // spin pos data initialized on device. [cx, cy, cz] for each.   [pd, r1,
    // r2, r2s, dB0]
    gPos_ = af::array(3, num, raw_loc_.data());
    raw_loc_.clear();

    // dB0
    gdB0s_ = af::array(1, num, raw_dB0_.data());
    gdB0s_ *= 267.5216941e3;  // into rad/s.
    raw_dB0_.clear();

    tx_sens_amp_ = af::array(nTxs_, num, raw_tx_B1_amp_.data());
    tx_sens_phi_ = af::array(nTxs_, num, raw_tx_B1_phase_.data());

    rx_sens_amp_ = af::array(nRxs_, num, raw_rx_B1_amp_.data());
    rx_sens_phi_ = af::array(nRxs_, num, raw_rx_B1_phase_.data());

    mx_ = af::constant(0, 1, num, f64);
    my_ = af::constant(0, 1, num, f64);
    mz_ = af::array(1, num, raw_rho_.data());

    r1_ = af::array(1, num, raw_r1_.data());
    r2_ = af::array(1, num, raw_r2_.data());
    pd_ = af::array(1, num, raw_rho_.data());

    rx_sens_cx_ = rx_sens_amp_ * af::cos(rx_sens_phi_);
    rx_sens_cy_ = rx_sens_amp_ * af::sin(rx_sens_phi_);
#endif
  }
}

// void engine::set_observer(const sp_cx_vec &rho) {
//  cout << "~~~~\n";
//#ifdef DENSE_MATRIX_COMPUTE
//  unified_spinsys_.det = rho.toDense(); // temporarily used.
//#else
//  unified_spinsys_.det = rho;
//#endif
//  levante_ernst_correction(unified_spinsys_.det);
//  cout << unified_spinsys_.det << "\n";
//}
//
// void engine::set_observer(const string expr) {
//  if (unified_spinsys_.p_sys != nullptr)
//    cout << "###\n";
//  else
//    cout << "@@@@\n";
//  set_observer(unified_spinsys_.p_sys->smart_state(expr));
//}

void engine::evolution(timeline dt, const seq_const &ctrl) {
  seq_step val;
  val.cur = ctrl;
  val.dt = timeline2ms(dt) * 1e-3;  // unit into s.
  seq_step_list_.push_back(val);

  double t = timeline2ms(dt) * 1e-3;  // unit into s.

  // OpenMP CPU case.
  if (physx_model_ == _quantum_cpu) {
#pragma omp parallel for
    for (int i = 0; i < (int)ensemble_.size(); i++)
      evolution_for_each(t, ctrl, ensemble_[i]);

    if (ctrl.acq.adc && ctrl.acq.last == true)
      raw_signal_.push_back(accu_signal());

    // if (ctrl.acq.adc)
    // raw_signal_.push_back(accu_signal());
  } else if (physx_model_ == _bloch_gpu) {
#ifdef ARRAYFIRE_USE
    af::array E1 = af::exp(-t * r1_);
    af::array E2 = af::exp(-t * r2_);

    // size_t id = (size_t)(ctrl.acq.index);
    int nRxs_ = 1;  // single coil.
    if (ctrl.acq.adc && ctrl.acq.index == 1) {
      af::array mx = af::tile(mx_, nRxs_, 1);
      af::array my = af::tile(my_, nRxs_, 1);

      af::array rx_fid_x = rx_sens_cx_ * mx - rx_sens_cy_ * my;
      af::array rx_fid_y = rx_sens_cy_ * mx + rx_sens_cx_ * my;

      af::array tmp_fid_x = af::sum(rx_fid_x, 1);
      af::array tmp_fid_y = af::sum(rx_fid_y, 1);

      // signal for each rx coil at this acq point.
      double *host_x = tmp_fid_x.host<double>();
      double *host_y = tmp_fid_y.host<double>();

      gpu_signal_transfer(host_x, host_y);

      delete[] host_x;
      delete[] host_y;
    }
     //if (ctrl.acq.adc) {
    //cout << ctrl.acq.index << "\n";
    //}

    af::array deltaB = gdB0s_;
    // if (norm(ctrl) == 0) {
    // allocate grad ctrl vector [x/y/z] on device.
    af::array grad(3, 1, ctrl.grad.v.data());
    // calculate total magnetic field variation due to gradient field, unit in
    // mT.
    af::array deltaGB = matmulTN(grad, gPos_);  // 1*nspins_
    deltaGB *= 267.5216941e3;                   // into rad/s.
    deltaB += deltaGB;
    //}

    af::array bufferMx, bufferMy, bufferMz;
    if (ctrl.rf_if) {
      af::array sens_amp = tx_sens_amp_;
      af::array sens_phi = tx_sens_phi_;

      // all coil have the same pulse pars.

      // only for 1H case (single spin).
      cd u(ctrl.rf[0].u[cx], ctrl.rf[0].u[cy]);
      sens_amp *= abs(u);
      sens_phi += phase_in_rad(u);
      deltaB += ctrl.rf[0].df * 2 * _pi;

      ux_ = sens_amp * af::cos(sens_phi);
      uy_ = sens_amp * af::sin(sens_phi);

      af::array sux = af::sum(ux_, 0);
      af::array suy = af::sum(uy_, 0);

      af::array rf_amp = af::sqrt(af::pow(sux, 2) + af::pow(suy, 2));
      af::array rf_phi = af::atan2(sux, suy);
      rf_phi *= -1;

      af::array Alpha = af::sqrt(af::pow(deltaB, 2) + af::pow(rf_amp, 2)) * t;
      af::array Beta = af::atan2(deltaB, rf_amp);

      af::array cosPhi = af::cos(rf_phi);
      af::array sinPhi = af::sin(rf_phi);
      af::array sinAlpha = af::sin(Alpha);
      af::array cosAlpha = af::cos(Alpha);
      af::array sinBeta = af::sin(Beta);
      af::array cosBeta = af::cos(Beta);
      af::array sinBeta2 = af::pow(sinBeta, 2);
      af::array cosBeta2 = af::pow(cosBeta, 2);

      af::array T1 = sinAlpha * sinPhi;
      af::array T2 = sinAlpha * cosPhi;
      af::array T3 = cosAlpha * sinPhi;
      af::array T4 = cosAlpha * cosPhi;
      af::array T5 = sinPhi * sinBeta;
      af::array T6 = cosPhi * sinBeta;
      af::array T7 = cosBeta * sinBeta;
      af::array T8 = cosAlpha * T7;
      af::array T9 = T4 + T1 * sinBeta;
      af::array T10 = cosBeta2 * cosPhi + sinBeta * (T1 + T4 * sinBeta);
      af::array T11 = cosBeta2 * sinPhi - sinBeta * (T2 - cosAlpha * T5);

      bufferMx = mx_ * (cosPhi * T9 + sinPhi * T11) -
                 my_ * (cosPhi * T11 - sinPhi * T9) +
                 mz_ * (cosBeta * (T2 - cosAlpha * T5) + cosBeta * T5);

      bufferMy = my_ * (cosPhi * T10 + sinPhi * (T3 - T2 * sinBeta)) -
                 mx_ * (sinPhi * T10 - cosPhi * (T3 - T2 * sinBeta)) +
                 mz_ * (cosBeta * (T1 + cosAlpha * T6) - cosBeta * T6);

      bufferMz = mx_ * (sinPhi * (T7 - T8) - cosBeta * T2) -
                 my_ * (cosPhi * (T7 - T8) + cosBeta * T1) +
                 mz_ * (sinBeta2 + cosAlpha * cosBeta2);

      mx_ = bufferMx * E2;
      my_ = bufferMy * E2;
      mz_ = bufferMz * E1 - (E1 - 1) * pd_;
    } else {
      af::array Alpha = deltaB * t;
      af::array sinAlpha = sin(Alpha);
      af::array cosAlpha = cos(Alpha);
      bufferMx = mx_ * cosAlpha - my_ * sinAlpha;
      bufferMy = my_ * cosAlpha + mx_ * sinAlpha;
      bufferMz = mz_;

      mx_ = bufferMx * E2;
      my_ = bufferMy * E2;
      mz_ = bufferMz * E1 - (E1 - 1) * pd_;
    }

    if (ctrl.acq.adc) {
      af::array mx = af::tile(mx_, nRxs_, 1);
      af::array my = af::tile(my_, nRxs_, 1);

      af::array rx_fid_x = rx_sens_cx_ * mx - rx_sens_cy_ * my;
      af::array rx_fid_y = rx_sens_cy_ * mx + rx_sens_cx_ * my;

      af::array tmp_fid_x = af::sum(rx_fid_x, 1);
      af::array tmp_fid_y = af::sum(rx_fid_y, 1);

      // signal for each rx coil at this acq point.
      double *host_x = tmp_fid_x.host<double>();
      double *host_y = tmp_fid_y.host<double>();

      gpu_signal_transfer(host_x, host_y);

      delete[] host_x;
      delete[] host_y;
    }

    if (ctrl.acq.adc && ctrl.acq.last == true)
      raw_signal_.push_back(accu_signal_bloch_gpu());
#endif
  }
}
void engine::gpu_signal_transfer(double *fidx, double *fidy) {
  Eigen::Map<vec> sx(fidx, nRxs_);
  Eigen::Map<vec> sy(fidy, nRxs_);
  cx_vec s(nRxs_);
  s.real() = sx;
  s.imag() = sy;
  raw_signal_tmp_.push_back(s);
}
void engine::evolution_for_each(double dt, const seq_const &ctrl,
                                each_spinsys &each) {
  each.L = each.L0 + each.R;

  // delay evolution.
  if (ctrl.delay_if) {
    // each.L += each.L0;
    each.rho = ssl::spinsys::step(each.rho, each.L, dt);
    // each.rho = ssl::spinsys::expmv(each.rho, -ci * each.L, dt, mat(1, 1), 1,
    // false);
    return;
  }

  // gradient part.
  double dw = each.dB;                               // mT.
  dw += (each.pos.cwiseProduct(ctrl.grad.v)).sum();  // mT.

  if (fabs(dw) > 1e-8) {
    dw *= 267.5216941e3;  // into rad/s.
    each.L += dw * each.Lz0;
  }

  // acquisition evolution (1st point).
  if (ctrl.acq.adc && ctrl.acq.index == 1) {
    each.sig = cx_vec(ctrl.acq.nps);
    each.sig[0] = projection(each.rho, unified_spinsys_.det);
    //    if (ctrl.acq.nps == 1) {
    //      if (ctrl.acq.last)
    //        each.rho = unified_spinsys_.rho0.toDense();
    //      return;
    //    }
  }

  // rf hamiltonians if any.
  if (ctrl.rf_if)
    for (size_t i = 0; i < ctrl.rf.size(); i++) {
      int ch = unified_spinsys_.rf_ctrl.channel_index(ctrl.rf[i].channel);
#ifdef DENSE_MATRIX_COMPUTE
      each.L += ctrl.rf[i].u[cx] * unified_spinsys_.rf_ctrl.Lx_dense[ch] +
                ctrl.rf[i].u[cy] * unified_spinsys_.rf_ctrl.Ly_dense[ch];
#else
      each.L += ctrl.rf[i].u[cx] * unified_spinsys_.rf_ctrl.Lx[ch] +
                ctrl.rf[i].u[cy] * unified_spinsys_.rf_ctrl.Ly[ch];
#endif
      double df = ctrl.rf[i].df;
      // cout << df << "####\n";
      if (df != 0) each.L += df * 2 * _pi * each.Lz0;
    }

  each.rho = ssl::spinsys::step(each.rho, each.L, dt);
  // each.rho = ssl::spinsys::expmv(each.rho, -ci * each.L, dt, mat(1, 1), 1,
  // false);

  // acquisition evolution (other points, acquire after the step evolution).
  if (ctrl.acq.adc) {
    each.sig[ctrl.acq.index] = projection(each.rho, unified_spinsys_.det);
    // cout << ctrl.acq.index << "###\n";
  }

  if (ctrl.acq.adc && ctrl.acq.last) {
    // each.rho = unified_spinsys_.rho0.toDense();
    // cout<< each.rho<<"\n\n";
  }
}

cx_vec engine::accu_signal_bloch_gpu() {
  int np = raw_signal_tmp_.size();
  int nrow = raw_signal_tmp_[0].size();
  cx_mat m(nrow, np);

#pragma omp parallel for
  for (int i = 0; i < np; i++) m.col(i) = raw_signal_tmp_[i];

  raw_signal_tmp_.clear();
  return m.row(0).transpose(); // currently only return 1st rx channel data.
}

cx_vec engine::accu_signal() {
  int npts = ensemble_[0].sig.size();
  cx_vec sig = cx_vec::Zero(npts);
  cx_vec ones = cx_vec::Ones(npts);

  vector<cx_vec> omp_sig(omp_core_num, sig);
#pragma omp parallel for
  for (int i = 0; i < (int)ensemble_.size(); i++) {
    int id = omp_get_thread_num();
    ensemble_[i].sig -= ones;
    omp_sig[id] += ensemble_[i].sig * ensemble_[i].pd;
  }

  for (int id = 0; id < omp_core_num; id++) sig += omp_sig[id];

  if (g_seq_param->acq_phase != 0) {
    // deal with the phase_acq
    for (int i = 0; i < sig.size(); i++) {
      cd a = fid_xy2amp(sig[i]);
      sig[i] = fid_amp2xy(cd(a.real(), a.imag() + g_seq_param->acq_phase));
    }
  }
  // optional
  //apodization(sig, 10);

  return sig;
  /*sol::table lines = g_lua->create_table();
  lines.add(sig.real());
  lines.add(sig.imag());
  plot("fid", line_series(lines));*/

  /*line fx(sig.real());
  line fy(sig.imag());

  (*g_lua)["_fid_x"] = fx;
  (*g_lua)["_fid_y"] = fy;
  g_lua->script("plot('title[fid]', _fid_x, _fid_y)");

  cx_vec spec = fft_1d(sig);
  line sx(spec.real());
  line sy(spec.imag());

  (*g_lua)["_spec_x"] = sx;
  (*g_lua)["_spec_y"] = sy;
  g_lua->script("plot('title[spec]', _spec_x, _spec_y)");*/
}
sol::object engine::process_signal() {
  ensemble_.clear(); // release memory for lots of spins.
  ssl_color_text("info", "proccesing acquisition data ...\n");

  int cols = raw_signal_[0].size();
  int rows = raw_signal_.size();
  ssl_color_text("info", "raw data size: " + to_string(rows) + "*" +
                             to_string(cols) + "\n");

  // lua table for raw data (fid, spec).
  sol::table t_raw = g_lua->create_table();

  // lua table for signal.
  sol::table t_fid = g_lua->create_table();
  sol::table t_fid_re = g_lua->create_table();
  sol::table t_fid_im = g_lua->create_table();
  sol::table t_fid_abs = g_lua->create_table();
  // lua table for spec.
  sol::table t_spec = g_lua->create_table();
  sol::table t_spec_re = g_lua->create_table();
  sol::table t_spec_im = g_lua->create_table();
  sol::table t_spec_abs = g_lua->create_table();

  cx_mat FID(rows, cols);
  cx_mat SPEC1D(rows, cols);

  for (size_t i = 0; i < raw_signal_.size(); i++) {
    cx_vec fid = raw_signal_[i];
    vec fid_re = fid.real();
    vec fid_im = fid.imag();
    vec fid_abs = fid.cwiseAbs();

    t_fid.add(fid);
    t_fid_re.add(fid_re);
    t_fid_im.add(fid_im);
    t_fid_abs.add(fid_abs);

    cx_vec spec = fft_1d(fid);
    vec spec_re = spec.real();
    vec spec_im = spec.imag();
    vec spec_abs = spec.cwiseAbs();
    t_spec.add(spec);
    t_spec_re.add(spec_re);
    t_spec_im.add(spec_im);
    t_spec_abs.add(spec_abs);

    FID.row(i) = fid.transpose();
    SPEC1D.row(i) = spec.transpose();
  }

  t_raw.set("fid", t_fid);
  t_raw.set("fid:re", t_fid_re);
  t_raw.set("fid:im", t_fid_im);
  t_raw.set("fid:abs", t_fid_abs);

  t_raw.set("spec", t_spec);
  t_raw.set("spec:re", t_spec_re);
  t_raw.set("spec:im", t_spec_im);
  t_raw.set("spec:abs", t_spec_abs);

  // for image.
  string time_s = sys_time();
  string folder = "raw_data_" + time_s;
  g_lua->script("os.execute('mkdir " + folder + "')");
  H5File file(folder + "/raw.h5", H5F_ACC_TRUNC);
  Group group1(file.createGroup("FID"));
  // raw data files.
  mat raw_re = FID.real();
  mat raw_im = FID.imag();
  mat raw_abs = FID.cwiseAbs();

  t_raw.set("FID", FID);
  t_raw.set("FID:re", raw_re);
  t_raw.set("FID:im", raw_im);
  t_raw.set("FID:abs", raw_abs);

  h5write(file, &group1, "FID:re", raw_re);
  h5write(file, &group1, "FID:im", raw_im);
  h5write(file, &group1, "FID:abs", raw_abs);

  //  (*g_lua)["_raw_fid_re"] = &raw_re;
  //  (*g_lua)["_raw_fid_im"] = &raw_im;
  //  (*g_lua)["_raw_fid_abs"] = &raw_abs;
  //
  //  g_lua->script("write('raw_fid_re.txt', _raw_fid_re)");
  //  g_lua->script("write('raw_fid_im.txt', _raw_fid_im)");
  //  g_lua->script("write('raw_fid_abs.txt', _raw_fid_abs)");

  // process the raw data if needed.
  cx_mat FID_post;
  if (g_phase_cycle_steps > 1) {
    size_t nrows = FID.rows();
    if (nrows % g_phase_cycle_steps != 0)
      throw std::runtime_error(
          "number of scans must be integral multiple of number of steps in the "
          "phase cycle!");

    FID_post = cx_mat(nrows / g_phase_cycle_steps, FID.cols());
    Eigen::RowVectorXcd sum_fid(FID.cols());
    for (size_t i = 0; i < FID_post.rows(); i++) {
      sum_fid.setZero();
      for (size_t j = 0; j < g_phase_cycle_steps; j++)
        sum_fid += FID.row(i * g_phase_cycle_steps + j);
      FID_post.row(i) = sum_fid;
    }
  } else
    FID_post = FID;

  cx_mat img = fft_2d(FID_post);
  Group group2(file.createGroup("IMG"));
  // raw data files.
  mat img_re = img.real();
  mat img_im = img.imag();
  mat img_abs = img.cwiseAbs();

  t_raw.set("IMG", img);
  t_raw.set("IMG:re", img_re);
  t_raw.set("IMG:im", img_im);
  t_raw.set("IMG:abs", img_abs);

  // h5write(file, &group2, "IMG", img);
  h5write(file, &group2, "IMG:re", img_re);
  h5write(file, &group2, "IMG:im", img_im);
  h5write(file, &group2, "IMG:abs", img_abs);
  h5write(file, &group2, "phase cycle steps", boost::lexical_cast<string>(g_phase_cycle_steps));

  Group group3(file.createGroup("SPEC"));
  // raw data files.
  mat spec_re = SPEC1D.real();
  mat spec_im = SPEC1D.imag();
  mat spec_abs = SPEC1D.cwiseAbs();

  t_raw.set("SPEC", SPEC1D);
  t_raw.set("SPEC:re", spec_re);
  t_raw.set("SPEC:im", spec_im);
  t_raw.set("SPEC:abs", spec_abs);

  h5write(file, &group3, "SPEC:re", spec_re);
  h5write(file, &group3, "SPEC:im", spec_im);
  h5write(file, &group3, "SPEC:abs", spec_abs);

  file.close();

  return t_raw;
}

}  // namespace physx
}  // namespace ssl
