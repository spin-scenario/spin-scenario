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
#include "seq_block.h"

namespace ssl {
namespace seq {

enum envelope_style {
  _amp_phase = 0,
  _ux_uy  // used in kernel calculations.
};
// RF characterization.
// RF carrier freq is tipically set equal to Larmor freq f(B0), +/- a specified freq offset delta_f.
// Here we use RF reference frame, that is using Larmor freq f(B0) minus the laboratory frequency of RF.
// for easy calculation.
struct RFChannel {
  cx_vec envelope;  // a slowly varying B1(t) or A(t).
  vec
      carrier;  // unit in rad/s. the specified freq offset delta_f. ref to p28 "handbook", delta-omega = gamma*B0 - omega_rf.
  string channel; // the isotope channel which is excited.
};

enum rf_loop_type {
  _var_freq = 0,
  _var_gain,
  _var_phase
};

class rf_pulse : public seq_block {
 public:
  rf_pulse();
  virtual ~rf_pulse();
  virtual rf_pulse *Clone() const = 0;

  virtual void get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const;
  virtual void evolution(int index = -1);
  virtual void write(ostream &ostr = cout) const;
  virtual void h5write(H5File &file, string abbr = "") const;
  virtual void plot() const;

  virtual mat get_shape() const;
  virtual string get_header() const;
  virtual void switch_rf_mode(string mode);
  virtual int switch2loop(int index);
  virtual vector<string> get_channels_str() const;
  void convert2(envelope_style mode);
  //virtual void set_modulated_gain(double g);
  virtual void set_gain(double g);
  inline envelope_style mode() const {
    return mode_;
  }

  inline double get_dt() const {
    return dt_;
  }
  inline size_t get_steps() const {
    return nsteps_;
  }
  inline size_t get_channels() const {
    return raw_data_.size();
  }
  inline size_t get_dims() const {
    return 2 * nsteps_ * raw_data_.size();
  }

  // ONLY USED FOR GRAPE OPTIMIZATION.
  vector<double> clone_raw_data() const;
  void update_raw_data(const double *vals);

  vector<cx_vec> export_signal() const;
  double sampling_freq() const;

  // rf power on all channels.  \int_0^tp {u^2(t)dt}
  virtual double rf_power() const;
 protected:
  virtual void assign();
 protected:
  size_t nsteps_;  ///< total segments of the rf pulse.
  size_t channels_;  ///< rf pulse channels, in MRI there is usually only 1H channel excitation.
  double dt_;  // time duration of each step. !!!!! unit in us.
  timeline tl_dt_;
  vector<RFChannel> raw_data_;  ///< <amp, phase> in <rad, rad> or <ux,uy> in <rad,rad>
  vector<RFChannel> raw_data0_;  ///< <amp, phase> in <rad, rad> or <ux,uy> in <rad,rad>

  envelope_style mode_;

  rf_loop_type loop_type_;

  double modulated_gain_; // amplitude gain.
  double modulated_phase_; // phase offset (rad).
  double modulated_freq_; // freq offset (Hz).

  vec modulated_freq_list_;
  vec modulated_gain_list_;
};
} /* namespace seq */
} /* namespace ssl */
