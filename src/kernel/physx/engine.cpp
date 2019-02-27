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

  string acc_str = "cpu";
  if (is_retrievable("acc", t)) {
    acc_str = retrieve_table_str("acc", t);
  }
  boost::to_lower(acc_str);
  if (acc_str == "gpu")
    physx_model_ = _quantum_gpu;
  else if (acc_str == "cpu") {
    physx_model_ = _quantum_cpu;
  }

  if (is_retrievable("spinsys", t)) {
    // general one spin system case.
    sol::object par = retrieve_table("spinsys", t);
    const spin_system &par_sys = par.as<const spin_system &>();
    unified_spinsys_.init(par_sys);
  } else {
    //g_lua->script("_sys = spin_system{B0 = '3 T', spins = '1H'}");
    g_lua->script("_sys = spin_system{spin = '1H'}");
    unified_spinsys_.init((*g_lua)["_sys"]);
  }

  p_phantom_ = nullptr;
  if (is_retrievable("phantom", t)) {
    string par = retrieve_table_str("phantom", t);
    p_phantom_ = new phantom(par.c_str());
  }

  init_ensemble(p_phantom_);
}

engine::~engine() {
}

void engine::init_ensemble(const phantom *p_phantom) {
  double ppm = 0;// 0.05;
  if (physx_model_ == _quantum_cpu) {
    if (!p_phantom) {
      ssl_color_text("warn", "no phantom specified, ONLY one spin system used.\n");
      each_spinsys each;
      ensemble_ = vector<each_spinsys>(1, each);
#ifdef DENSE_MATRIX_COMPUTE
      ensemble_[0].pos = vec3::Zero();
      ensemble_[0].rho = unified_spinsys_.rho0.toDense();
      ensemble_[0].L0 =
          unified_spinsys_.L0.toDense(); //////////////////////////////////////////////////////////////////////////
      ensemble_[0].Lz0 = unified_spinsys_.Lz0.toDense();
      ensemble_[0].R = unified_spinsys_.R.toDense();
      ensemble_[0].pd = 1;
#else
      ensemble_[0].pos = vec3::Zero();
      ensemble_[0].rho = unified_spinsys_.rho0;
      ensemble_[0].L0 = unified_spinsys_.L0; //////////////////////////////////////////////////////////////////////////
      ensemble_[0].Lz0 = unified_spinsys_.Lz0;
      ensemble_[0].R = unified_spinsys_.R;
      ensemble_[0].pd = 1;
#endif
      return;
    }

    int num = p_phantom->isochromats_.size();
    each_spinsys each;
    ensemble_ = vector<each_spinsys>(num, each);
//#pragma omp parallel for
    for (int i = 0; i < num; i++) {
      //int id = omp_get_thread_num();
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
      ensemble_[i].dB = p_phantom->isochromats_[i].data[dB0] * ppm * unified_spinsys_.B0 / 1e3; // into mT
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
      ensemble_[i].dB = p_phantom->isochromats_[i].data[dB0] * ppm * unified_spinsys_.B0 / 1e3; // into mT.
      ensemble_[i].pd = p_phantom->isochromats_[i].data[pd];
#endif
    }
  }
  if (physx_model_ == _quantum_gpu) {
#ifdef ARRAYFIRE_COMPUTE
    af::setDevice(0);
    af::info();

    int n =  p_phantom->isochromats_.size();
    int dim = unified_spinsys_.rho0.size();

    vector<cd> raw_rhos;
    cx_vec rho0 = unified_spinsys_.rho0.toDense();
    vector<cd> raw_rho0 = vector<cd>(&rho0.data()[0], &rho0.data()[dim]);

    vector<double> raw_pos;

    vector<cd> raw_Lz0s;
    cx_mat Lz0 = unified_spinsys_.Lz0.toDense();
    vector<cd> raw_Lz0 = vector<cd>(&Lz0.data()[0], &Lz0.data()[dim*dim]);

    vector<cd> raw_L0s;
    cx_mat L0 = unified_spinsys_.L0.toDense();
    vector<cd> raw_L0 = vector<cd>(&L0.data()[0], &L0.data()[dim*dim]);

    vector<cd> raw_Rs;

    vector<double> raw_dB0s;

    for (int i = 0; i < n; i++) {
        // initial state.
        std::copy(raw_rho0.begin(), raw_rho0.end(), std::back_inserter(raw_rhos));

        // position.
        vec3 p = p_phantom->isochromats_[i].position();
        vector<double> pv(&p.data()[0], &p.data()[3]);
        std::copy(pv.begin(), pv.end(), std::back_inserter(raw_pos));

        //std::copy(raw_Lz0.begin(), raw_Lz0.end(), std::back_inserter(raw_Lz0s));

        std::copy(raw_L0.begin(), raw_L0.end(), std::back_inserter(raw_L0s));

        cx_mat R(dim, dim);
        R.setZero();
        R(2, 0) = ci*p_phantom->isochromats_[i].data[r1];
        R(1, 1) = -ci*p_phantom->isochromats_[i].data[r2];
        R(2, 2) = -ci*p_phantom->isochromats_[i].data[r1];
        R(3, 3) = -ci*p_phantom->isochromats_[i].data[r2];
        vector<cd> raw_R = vector<cd>(&R.data()[0], &R.data()[dim*dim]);
        std::copy(raw_R.begin(), raw_R.end(), std::back_inserter(raw_Rs));

        double dB = p_phantom->isochromats_[i].data[dB0] * ppm * unified_spinsys_.B0 / 1e3; // into mT.
        raw_dB0s.push_back(dB);
    }

    af_ensemble_.n = n;
    af_ensemble_.rho = af::array(dim, n, (af::cdouble*)raw_rhos.data());
    af_ensemble_.pos = af::array(3, n, (double*)raw_pos.data());
    af_ensemble_.L0 = af::array(dim, dim, n, (af::cdouble*)raw_L0s.data()); // to be only 1.
    af_ensemble_.Lz0 = af::array(dim, dim, (af::cdouble*)raw_Lz0.data()); // to be only 1.
    af_ensemble_.L = af::constant(0, af::dim4(dim, dim, n), c64);
    af_ensemble_.R = af::array(dim, dim, n, (af::cdouble*)raw_Rs.data());

    af_ensemble_.dB0 = af::array(1, n, (double*)raw_dB0s.data());

#ifdef DENSE_MATRIX_COMPUTE
    cx_vec det = unified_spinsys_.det;
#else
    cx_vec det = unified_spinsys_.det.toDense();
#endif
    vector<cd> raw_det = vector<cd>(&det.data()[0], &det.data()[dim]);
    af_ensemble_.det = af::array(dim, 1, (af::cdouble*)raw_det.data());


    int channels = unified_spinsys_.rf_ctrl.channels;
    vector<cd> raw_Lxs, raw_Lys;
    for (int i = 0; i < channels; i++) {
        cx_mat& Lx = unified_spinsys_.rf_ctrl.Lx_dense[i];
        vector<cd> raw_Lx = vector<cd>(&Lx.data()[0], &Lx.data()[dim*dim]);
        cx_mat& Ly = unified_spinsys_.rf_ctrl.Ly_dense[i];
        vector<cd> raw_Ly = vector<cd>(&Ly.data()[0], &Ly.data()[dim*dim]);
        std::copy(raw_Lx.begin(), raw_Lx.end(), std::back_inserter(raw_Lxs));
        std::copy(raw_Ly.begin(), raw_Ly.end(), std::back_inserter(raw_Lys));
    }

    af_ensemble_.Lx = af::array(dim, dim, channels, (af::cdouble*)raw_Lxs.data());
    af_ensemble_.Ly = af::array(dim, dim, channels, (af::cdouble*)raw_Lys.data());

    //af_ensemble_.dB0 *= 1e5;
    /*af_print(af_ensemble_.det);
    af_print(af_ensemble_.dB0);
    af_print(af_ensemble_.rho);
    af_print(af_ensemble_.pos);
    af_print(af_ensemble_.Lz0);
    af_print(af_ensemble_.R);
    af_print(af_ensemble_.L);
    af_print(af_ensemble_.Lx);
    af_print(af_ensemble_.Ly);
    exit(0);*/
#endif
  }
}

//void engine::set_observer(const sp_cx_vec &rho) {
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
//void engine::set_observer(const string expr) {
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
    omp_set_num_threads(omp_core_num);
#pragma omp parallel for
    for (int i = 0; i < (int) ensemble_.size(); i++)
      evolution_for_each(t, ctrl, ensemble_[i]);

    if (ctrl.acq.adc && ctrl.acq.last == true)
      raw_signal_.push_back(accu_signal());

    //if (ctrl.acq.adc)
    //raw_signal_.push_back(accu_signal());
  }

  // ArrayFire GPU case.
  if (physx_model_ == _quantum_gpu) {
#ifdef ARRAYFIRE_COMPUTE
    af_ensemble_.L = af_ensemble_.R;
    // delay evolution.
    if (ctrl.delay_if) {
        af_ensemble_.L += af_ensemble_.L0;
        af_step(af_ensemble_.n, af_ensemble_.L, af_ensemble_.rho, t);
        //each.rho = ssl::spinsys::step(each.rho, each.L, dt);
        return;
    }

    // gradient part.
    af::array dw = af_ensemble_.dB0; // mT.
    if (ctrl.grad.v.norm()) {
        af::array grad(3, 1, (double*)ctrl.grad.v.data()); // allocate grad ctrl vector [x/y/z] on device.
        dw += af::matmulTN(grad, af_ensemble_.pos); // 1*nspins_
    }

    dw *= 267.5216941e3; // into rad/s.

    gfor(af::seq k, af_ensemble_.n) {
        af_ensemble_.L(af::span, af::span, k) += dw(af::span, k)*af_ensemble_.Lz0;
    }

    // acquisition evolution (1st point).
    if (ctrl.acq.adc && ctrl.acq.index == 1) {
        af_ensemble_.sig = af::constant(0, af::dim4(ctrl.acq.nps, af_ensemble_.n), c64);
        af_ensemble_.sig(0, af::span) = af::matmulTN(af_ensemble_.det, af_ensemble_.rho);
    }

    // rf hamiltonians if any.
    if (ctrl.rf_if)
        for (int i = 0; i < ctrl.rf.size(); i++) {
            int ch = unified_spinsys_.rf_ctrl.channel_index(ctrl.rf[i].channel);
            gfor(af::seq k, af_ensemble_.n) {
                af_ensemble_.L(af::span, af::span, k) += ctrl.rf[i].u[cx] * af_ensemble_.Lx(af::span, af::span, ch) + ctrl.rf[i].u[cy] * af_ensemble_.Ly(af::span, af::span, ch);
            }
        }

    af_step(af_ensemble_.n, af_ensemble_.L, af_ensemble_.rho, t);

    if (ctrl.acq.adc)
        af_ensemble_.sig(ctrl.acq.index, af::span) = af::matmulTN(af_ensemble_.det, af_ensemble_.rho);

    if (ctrl.acq.adc && ctrl.acq.last == true) {
        af::array fid = af::sum(af_ensemble_.sig, 1); // sum along rows.
        af::cdouble* host_fid = fid.host<af::cdouble>();
        cx_vec v(ctrl.acq.nps);
        for (int i = 0; i < ctrl.acq.nps; i++)
            v[i] = cd(host_fid[i].real, host_fid[i].imag);
        raw_signal_.push_back(v);
    }
#endif
  }
}

void engine::evolution_for_each(double dt, const seq_const &ctrl, each_spinsys &each) {
  each.L = each.L0 + each.R;

  // delay evolution.
  if (ctrl.delay_if) {
    //each.L += each.L0;
    each.rho = ssl::spinsys::step(each.rho, each.L, dt);
    //each.rho = ssl::spinsys::expmv(each.rho, -ci * each.L, dt, mat(1, 1), 1, false);
    return;
  }

  // gradient part.
  double dw = each.dB; // mT.
  dw += (each.pos.cwiseProduct(ctrl.grad.v)).sum(); // mT.

  if (fabs(dw) > 1e-8) {
    dw *= 267.5216941e3; // into rad/s.
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
      each.L += ctrl.rf[i].u[cx] * unified_spinsys_.rf_ctrl.Lx_dense[ch]
          + ctrl.rf[i].u[cy] * unified_spinsys_.rf_ctrl.Ly_dense[ch];
#else
      each.L += ctrl.rf[i].u[cx] * unified_spinsys_.rf_ctrl.Lx[ch] + ctrl.rf[i].u[cy] * unified_spinsys_.rf_ctrl.Ly[ch];
#endif
      double df = ctrl.rf[i].df;
      //cout << df << "####\n";
      if (df != 0)
        each.L += df * 2 * _pi * each.Lz0;
    }

  each.rho = ssl::spinsys::step(each.rho, each.L, dt);
  //each.rho = ssl::spinsys::expmv(each.rho, -ci * each.L, dt, mat(1, 1), 1, false);

  // acquisition evolution (other points, acquire after the step evolution).
  if (ctrl.acq.adc) {
    each.sig[ctrl.acq.index] = projection(each.rho, unified_spinsys_.det);
    //cout << ctrl.acq.index << "###\n";
  }

  if (ctrl.acq.adc && ctrl.acq.last) {
    //each.rho = unified_spinsys_.rho0.toDense();
    //cout<< each.rho<<"\n\n";
  }
}

cx_vec engine::accu_signal() {
  int npts = ensemble_[0].sig.size();
  cx_vec sig = cx_vec::Zero(npts);
  cx_vec ones = cx_vec::Ones(npts);

  vector<cx_vec> omp_sig(omp_core_num, sig);
#pragma omp parallel for
  for (int i = 0; i < (int) ensemble_.size(); i++) {
    int id = omp_get_thread_num();
    ensemble_[i].sig -= ones;
    omp_sig[id] += ensemble_[i].sig * ensemble_[i].pd;
  }

  for (int id = 0; id < omp_core_num; id++)
    sig += omp_sig[id];

  if (g_seq_param->acq_phase != 0) {
    // deal with the phase_acq
    for (int i = 0; i < sig.size(); i++) {
      cd a = fid_xy2amp(sig[i]);
      sig[i] = fid_amp2xy(cd(a.real(), a.imag() + g_seq_param->acq_phase));
    }
  }
  // optional
  //apodization(sig, 5);

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
  ssl_color_text("info", "proccesing acquisition data ...\n");

  int cols = raw_signal_[0].size();
  int rows = raw_signal_.size();
  ssl_color_text("info", "raw data size: " + to_string(rows) + "*" + to_string(cols) + "\n");

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

//    double phi0=deg2rad(117);
//    double phi_i = 0.5*phi0*(double)(i*i+i+2);
//
//    for(size_t j=0;j<fid.size();j++)
//    {
//      cd tmp=xy2amp(fid[j]);
//      tmp=cd(tmp.real(), tmp.imag()-phi_i);
//      fid[i]=amp2xy(tmp);
//    }


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



  // for iamge.
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

  h5write(file, &group1, "fid_re", raw_re);
  h5write(file, &group1, "fid_im", raw_im);
  h5write(file, &group1, "fid_abs", raw_abs);

//  (*g_lua)["_raw_fid_re"] = &raw_re;
//  (*g_lua)["_raw_fid_im"] = &raw_im;
//  (*g_lua)["_raw_fid_abs"] = &raw_abs;
//
//  g_lua->script("write('raw_fid_re.txt', _raw_fid_re)");
//  g_lua->script("write('raw_fid_im.txt', _raw_fid_im)");
//  g_lua->script("write('raw_fid_abs.txt', _raw_fid_abs)");



  cx_mat img = fft_2d(FID);
  Group group2(file.createGroup("IMG"));
  // raw data files.
  mat img_re = img.real();
  mat img_im = img.imag();
  mat img_abs = img.cwiseAbs();

  t_raw.set("IMG", img);
  t_raw.set("IMG:re", img_re);
  t_raw.set("IMG:im", img_im);
  t_raw.set("IMG:abs", img_abs);

  //h5write(file, &group2, "IMG", img);
  h5write(file, &group2, "IMG:re", img_re);
  h5write(file, &group2, "IMG:im", img_im);
  h5write(file, &group2, "IMG:abs", img_abs);

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
  // 1d nmr.
//  cx_vec fid_apo = fid.row(0).transpose();
//  apodization(fid_apo, 5);
//
//  cx_vec fid_zero(8192);
//  fid_zero.setZero();
//  fid_zero.head(fid_apo.size()) = fid_apo;
//
//  cx_vec spec = fft_1d(fid_zero);
//
//  vec hz = vec::LinSpaced(spec.size(), -1e3 * g_seq_param->sw / 2, 1e3 * g_seq_param->sw / 2);
//  utility::line amp(hz, spec.cwiseAbs());
//  (*g_lua)["_amp"] = amp;
//  g_lua->script("plot('title#spec amp# gnuplot#set xrange [] reverse#', _amp)");


//  double max_noise = fid.cwiseAbs().maxCoeff()/100;
//  cx_mat noise(rows, cols);
//  noise.setRandom();
//  //fid+=noise;
//

//  utility::map map(img_abs.matrix());
//  (*g_lua)["_map"] = map;
//  g_lua->script("plot('title#image# gnuplot#set size ratio -1\\n set palette gray#', _map)");
//  (*g_lua)["_raw_img_re"] = &img_re;
//  (*g_lua)["_raw_img_im"] = &img_im;
//  (*g_lua)["_raw_img_abs"] = &img_abs;
//
//  g_lua->script("ssl.write('raw_img_re.txt', _raw_img_re)");
//  g_lua->script("ssl.write('raw_img_im.txt', _raw_img_im)");
//  g_lua->script("ssl.write('raw_img_abs.txt', _raw_img_abs)");

  // 1D-FFT array.
  /*cx_mat spec = fid;
  for (int i = 0; i < spec.rows(); i++)
      spec.row(i) = fft_1d(fid.row(i).transpose()).transpose();

  mat spec_abs = spec.cwiseAbs();

  ssl::utility::array amp_mat;
  amp_mat = 20 * (spec_abs.array() + 1e-6).log10();

  mat spec_phase = spec_abs;
  for (int i = 0; i < spec.rows(); i++)
      for (int j = 0; j < spec.cols(); j++)
          spec_phase(i, j) = phase_in_degree(spec(i, j));

  utility::map map1(spec_abs.matrix());
  (*g_lua)["_map1"] = map1;
  g_lua->script("plot('title<spec abs map> ylabel<Freq>', _map1)");

  utility::map map2(amp_mat.matrix());
  (*g_lua)["_map2"] = map2;
  g_lua->script("plot('title<spec 20log map> ylabel<Freq>', _map2)");

  utility::map map3(spec_phase.matrix());
  (*g_lua)["_map3"] = map3;
  g_lua->script("plot('title<spec phase map> ylabel<Freq>', _map3)");

  int id = spec.cols() / 2;
  utility::line amp(spec.col(id).cwiseAbs());
  (*g_lua)["_amp"] = amp;
  g_lua->script("plot('title<spec amp>', _amp)");*/


  /*vec m3 = fid.cwiseAbs().row(0).transpose();
  vec m1 = fid.real().row(0).transpose();
  vec m2 = fid.imag().row(0).transpose();

  sol::table lines = g_lua->create_table();
  lines.add(m1);
  lines.add(m2);
  lines.add(m3);
  plot("title<xy profile>", line_series(lines));*/



  //cx_vec sig = raw_signal_[0];
  //cout << sig[0] << "\n";
  //// plot.
  //sol::table lines = g_lua->create_table();
  //vec re = sig.real();
  //vec im = sig.imag();
  //vec abs = sig.cwiseAbs();

  //lines.add(re);
  //lines.add(im);
  //lines.add(abs);

  //string fig_spec;
  //fig_spec = "title<fid>";
  //plot(fig_spec, line_series(lines));

  //sol::table lines1 = g_lua->create_table();
  //cx_vec spec = fft_1d(sig);
  // re = spec.real();
  // im = spec.imag();
  // abs = spec.cwiseAbs();

  //lines1.add(re);
  //lines1.add(im);
  //lines1.add(abs);

  //string fig_spec1;
  //fig_spec1 = "title<spec>";
  //plot(fig_spec1, line_series(lines1));


//   double max_amp = fid.cwiseAbs().maxCoeff();
//   fid /= max_amp;

  //cx_mat img = fft_2d(fid);

  ////-----------------------PLOT-------------------------------
//   mat fid_re = fid.real();
//   mat fid_im = fid.imag();
  //mat fid_abs = fid.cwiseAbs();

  //mat img_re = img.real();
  //mat img_im = img.imag();
  //mat img_abs = img.cwiseAbs();

//   (*g_lua)["_raw_sig"] = &fid;
//   (*g_lua)["_raw_sig_re"] = &fid_re;
//   (*g_lua)["_raw_sig_im"] = &fid_im;
  //(*g_lua)["_raw_sig_abs"] = &fid_abs;

  //(*g_lua)["_raw_img"] = &img;
  //(*g_lua)["_raw_img_re"] = &img_re;
  //(*g_lua)["_raw_img_im"] = &img_im;
  //(*g_lua)["_raw_img_abs"] = &img_abs;

  //g_lua->script("os.execute('mkdir signal')");

//   g_lua->script("write('signal/raw_sig', _raw_sig)");
//   g_lua->script("write('signal/raw_sig_re', _raw_sig_re)");
//   g_lua->script("write('signal/raw_sig_im', _raw_sig_im)");
  //g_lua->script("write('signal/raw_sig_abs', _raw_sig_abs)");

  //g_lua->script("write('signal/raw_img', _raw_img)");
  //g_lua->script("write('signal/raw_img_re', _raw_img_re)");
  //g_lua->script("write('signal/raw_img_im', _raw_img_im)");
  //g_lua->script("write('signal/raw_img_abs', _raw_img_abs)");

//   string gnu_cmd = "gnuplot<set yrange [" + to_string(fid.rows()) + ":1]>";
//   (*g_lua)["_map_img"] = utility::map(img.cwiseAbs());
//   g_lua->script("plot('title<img> " + gnu_cmd + "', _map_img)");
//   (*g_lua)["_map_fid"] = utility::map(fid.cwiseAbs());
//   g_lua->script("plot('title<fid> " + gnu_cmd + "', _map_fid)");
}

}
}
