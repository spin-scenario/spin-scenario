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

#include "shaped_rf.h"
#include <unsupported/Eigen/Splines>
namespace ssl {
namespace seq {
class spline {
 public:
  spline(vec const &x_vec, vec const &y_vec)
      : x_min(x_vec.minCoeff()),
        x_max(x_vec.maxCoeff()),
        spline_(Eigen::SplineFitting<Eigen::Spline<double, 1>>::Interpolate(
            y_vec.transpose(),
            std::min<int>(x_vec.rows() - 1, 3),
            scaled_values(x_vec))) {}
  double operator()(double x) const {
    return spline_(scaled_value(x))(0);
  }
 private:
  double scaled_value(double x) const {
    return (x - x_min) / (x_max - x_min);
  }
  Eigen::RowVectorXd scaled_values(vec const &x_vec) const {
    return x_vec.unaryExpr([this](double x) {
      return scaled_value(x);
    }).transpose();
  }
  double x_min;
  double x_max;
  Eigen::Spline<double, 1> spline_;
  // http://stackoverflow.com/questions/29822041/eigen-spline-interpolation-how-to-get-spline-y-value-at-arbitray-point-x
};

shaped_rf::shaped_rf() {
}

shaped_rf::~shaped_rf() {
  // TODO Auto-generated destructor stub
}

void shaped_rf::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "shaped rf pulse"); // the config table has been initialized in seq_block().

  std::string s_name = name();

  double width = retrieve_config_table_double("width"); // required.
  timer_.width = ms2timeline(width);

  std::string channel = "1H"; // if not specified, use 1H as default.
  if (is_retrievable("channel")) {
    channel = retrieve_config_table_str("channel"); // required.
  }

  std::vector<std::string> par_vec;
  boost::split(par_vec, channel, boost::is_any_of(", |"), boost::token_compress_on);
  channels_ = par_vec.size();

  mode_ = _amp_phase;
  if (is_retrievable("mode")) {
    std::string mode = retrieve_config_table_str("mode");
    if (mode == "amp/phase")
      mode_ = _amp_phase;
    else if (mode == "ux/uy")
      mode_ = _ux_uy;
    else {
      std::string s = "unknown shaped pulse mode, either amp/phase or ux/uy!";
      throw std::runtime_error(s.c_str());
    }
  }

  std::string waveform = retrieve_config_table_str("pattern");

  std::map<std::string, rf_pattern>::const_iterator iter;
  iter = g_rf_pattern.find(waveform);
  if (iter != g_rf_pattern.end()) {
    if (is_retrievable("step"))  // external file do not need to assign this par.
      nsteps_ = (size_t) (retrieve_config_table_int("step"));

    timeline dt = ms2timeline(width) / nsteps_;
    dt_ = timeline2ms(dt) * 1e3;  // us.
    timer_.width = dt * nsteps_;
    tl_dt_ = dt;

    double max_amp = 1;
    if (is_retrievable("max_amp"))
      max_amp = retrieve_config_table("max_amp").as<double>();
    set_shape(par_vec, iter->second, max_amp);
  } else {
    set_shape(par_vec, waveform); //external file.
    timeline dt = ms2timeline(width) / nsteps_;
    dt_ = timeline2ms(dt) * 1e3;  // us.
    timer_.width = dt * nsteps_;
    tl_dt_ = dt;
  }

  // deal with scaling if required.
  if (is_retrievable("beta")) {
    double beta = retrieve_config_table_double("beta"); // unit in degree.
    beta = deg2rad(beta);// into rad.

    convert2(_amp_phase);

    for (size_t k = 0; k < raw_data_.size(); k++) {

      vec amp = raw_data_[k].envelope.real();

      raw_data_[k].envelope.real() *= beta / (amp.sum() * dt_ * 1e-6);

      //double ratio = beta/(_pi*2)/(amp.sum() * dt_ * 1e-6);
      //raw_data_[k].envelope.real() *= ratio;
    }

  }

  if (timer_.width != ms2timeline(width)) {
    std::string s = str(boost::format("%s %s %s ms.\n") % name() % "rf duration auto adjusted to " % width_in_ms());
    ssl_color_text("warn", s);
  }
  timer_.keys = tlvec::LinSpaced(nsteps_ + 1, 0, timer_.width);

  s_name += " - " + waveform;

  config_table_.set("name", s_name);
  rf_pulse::assign();
}

void shaped_rf::set_shape(std::vector<std::string> channels, const mat &m) {
  raw_data_.clear(); // NOTE!!!
  nsteps_ = m.rows(); //!!! auto assign step from files.
  if (m.cols() != (int) (2 * channels_)) {
    std::string s = "shaped pulse data does not match the specified channel!";
    throw std::runtime_error(s.c_str());
  }
  for (size_t i = 0; i < channels_; i++) {
    RFChannel chl;
    chl.envelope = cx_vec::Zero(nsteps_);
    chl.carrier = vec::Zero(nsteps_);
    chl.envelope.real() = 2 * _pi * m.col(i * 2); // into rad.
    chl.envelope.imag() = m.col(i * 2 + 1);
    if (mode_ == _ux_uy)
      chl.envelope.imag() *= 2 * _pi; // uy in rad.
    if (mode_ == _amp_phase)
      chl.envelope.imag() *= _pi / 180; // phase in rad.

    chl.channel = channels[i];
    raw_data_.push_back(chl);
  }
  raw_data0_ = raw_data_;
}

void shaped_rf::set_shape(std::vector<std::string> channels, std::string file) {
  raw_data_.clear(); // NOTE!!!
  mat data;
  if (boost::ends_with(file.c_str(), ".h5")) {
    hid_t file_id;
    herr_t status;
    file_id = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    g_h5_string.clear();
    status = H5Literate(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func, NULL);
    status = H5Fclose(file_id);
    if (status < 0)
      return;

    H5File h5file;
    h5file.openFile(file, H5F_ACC_RDWR);
    data = h5read_mat(h5file, "/" + g_h5_string[0] + "/shape");
    h5file.close();
  } else if (boost::ends_with(file.c_str(), ".RF"))
    data = eigen_read(file);
  else // YACAS expr.
  {
    nsteps_ = (size_t) (retrieve_config_table_int("step")); // required.
    std::string pattern_expr = yacas_evaluate(file); // in case of D(t) style.

    data = mat(nsteps_, 2 * channels.size());
    vec time = vec::LinSpaced(nsteps_, 0, width_in_ms());
    yacas_evaluate("f(t):=" + pattern_expr);
    for (int i = 0; i < nsteps_; i++) {
      data(i, 0) = std::stod(yacas_evaluate("N(f(" + std::to_string(time(i)) + "))"));
      data(i, 1) = 0; // phase.
    }
    // TODO: to be added for multi-expr of diff channels.
  }

  if (data.size() == 0) {
    std::string s = "failed to read shaped pulse data from file: " + file;
    throw std::runtime_error(s.c_str());
  }
  set_shape(channels, data);
}

void shaped_rf::set_shape(std::vector<std::string> channels, rf_pattern p, double maxamp) {
  raw_data_.clear(); // NOTE!!!
  size_t steps = nsteps_;
  maxamp *= 2 * _pi; // into rad.
  for (size_t i = 0; i < channels_; i++) {
    RFChannel chl;
    chl.envelope = cx_vec::Zero(steps);
    chl.carrier = vec::Zero(steps);
    double width = width_in_ms();
    switch (p) {
      case _sinc: {
        double lobes = 5;
        if (is_retrievable("lobe")) {
          lobes = retrieve_config_table_double("lobe"); // required.
        }
        lobes = lobes + 1;
        double t0 = (double) width / lobes;  // ms
        vec vt = vec::LinSpaced(steps, -width / 2, width / 2);
        for (size_t i = 0; i < steps; i++) {
          double amp = t0 * sin(_pi * vt(i) / t0) / (_pi * vt(i));
          chl.envelope[i] =
              cd(fabs(amp) * maxamp, amp > 0 ? 0 : deg2rad(180));  // <amplitude (Hz), phase (deg)>, x-pulse.
        }
      }
        break;
      case _rect: {
        vec rect = window_function(utility::_rect, steps);
        chl.envelope.real() = rect * maxamp;
      }
        break;
      case _gauss: {
        vec gau = window_function(utility::_gauss, steps);
        for (size_t i = 0; i < steps; i++) {
          double amp = gau[i];
          chl.envelope[i] =
              cd(fabs(amp) * maxamp, amp > 0 ? 0 : deg2rad(180));  // <amplitude (Hz), phase (deg)>, x-pulse.
        }
      }
        break;
      case _hamming: {
        vec ham = window_function(utility::_hamming, steps);
        for (size_t i = 0; i < steps; i++) {
          double amp = ham[i];
          chl.envelope[i] =
              cd(fabs(amp) * maxamp, amp > 0 ? 0 : deg2rad(180));  // <amplitude (Hz), phase (deg)>, x-pulse.
        }
      }
        break;
      case _rand: {
        chl.envelope.real() = vec::Random(steps).cwiseAbs() * maxamp;
        chl.envelope.imag() = vec::Random(steps) * _pi;
      }
        break;
      case _rand_spline: {
        vec amp = vec::Random(steps).cwiseAbs() * maxamp;
        vec phase = vec::Random(steps).cwiseAbs() * 360;
        size_t span = 15;
        size_t N = (steps - 1) / span + 2;
        vec xvals(N), amp_vals(N), phase_vals(N);
        size_t i = 0, j = 0;
        for (i = 0; i < steps; i += span) {
          xvals[j] = i;
          amp_vals[j] = amp[i];
          phase_vals[j] = phase[i];
          j++;
        }
        if (i != steps - 1) {
          xvals[j] = steps - 1;
          amp_vals[j] = amp[steps - 1];
          phase_vals[j] = phase[steps - 1];
        }

        spline s1(xvals, amp_vals);
        spline s2(xvals, phase_vals);
        for (i = 0; i < steps; i++) {
          amp[i] = s1(i);
          phase[i] = s2(i);
          if (amp[i] < 0)
            amp[i] = 0;
          if (amp[i] > maxamp)
            amp[i] = maxamp;
          if (phase[i] > 360)
            phase[i] = 360;
          if (phase[i] < 0)
            phase[i] = 0;
        }
        chl.envelope.real() = amp;
        chl.envelope.imag() = phase;
        chl.envelope.imag() *= _pi / 180; // phase into rad.
      }
        break;
      case _user_defined: {

      }
        break;
      default:break;
    }
    chl.channel = channels[i];
    raw_data_.push_back(chl);
  }

  if (mode_ == _ux_uy) {
    mode_ = _amp_phase; // by default in amp/phase mode.
    convert2(_ux_uy);
  }
}

std::map<std::string, rf_pattern> rf_pattern_map() {
  std::map<std::string, rf_pattern> map;
  map.insert(std::pair<std::string, rf_pattern>("sinc", _sinc));
  map.insert(std::pair<std::string, rf_pattern>("gauss", _gauss));
  map.insert(std::pair<std::string, rf_pattern>("hamming", _hamming));
  map.insert(std::pair<std::string, rf_pattern>("rand", _rand));
  map.insert(std::pair<std::string, rf_pattern>("rect", _rect));
  map.insert(std::pair<std::string, rf_pattern>("rand_spline", _rand_spline));
  map.insert(std::pair<std::string, rf_pattern>("user_defined", _user_defined));
  return map;
}

} /* namespace seq */
} /* namespace ssl */
