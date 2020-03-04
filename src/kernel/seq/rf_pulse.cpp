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

#include "rf_pulse.h"
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;
#include <chrono>
using namespace chrono;
namespace ssl {
namespace seq {

rf_pulse::rf_pulse()
    : nsteps_(0),
      channels_(1),
      dt_(0),
      mode_(_amp_phase) {
  modulated_gain_ = 1;
  modulated_freq_ = 0;
}

rf_pulse::~rf_pulse() {
  // TODO Auto-generated destructor stub
}
void rf_pulse::assign() {
  config_table_.set("category", _rf);

  if (is_retrievable("gain")) {
    sol::object obj = retrieve_config_table("gain");
    if (obj.get_type() == sol::type::number)
      modulated_gain_ = obj.as<double>();

    if (obj.get_type() == sol::type::string) {
      colon_sep val;
      if (parse(obj.as<string>(), val)) {
        modulated_gain_list_ = vec::LinSpaced(val.num, val.a, val.b);
        loop_ctrl_.loop_count = modulated_gain_list_.size();
        loop_type_ = _var_gain;
      }
    }
  }

  if (is_retrievable("freq")) {
    sol::object obj = retrieve_config_table("freq");
    if (obj.get_type() == sol::type::number)
      modulated_freq_ = obj.as<double>();

    if (obj.get_type() == sol::type::string) {
      colon_sep val;
      if (parse(obj.as<string>(), val)) {
        modulated_freq_list_ = vec::LinSpaced(val.num, val.a, val.b);
        loop_ctrl_.loop_count = modulated_freq_list_.size();
        loop_type_ = _var_freq;
      }
    }
  }

  if (is_retrievable("phase")) {
    sol::object obj = retrieve_config_table("phase");
    if (obj.get_type() == sol::type::number) {
      loop_ctrl_.loop_count = obj.as<int>();
      loop_type_ = _var_phase;
    }
  }

  seq_block::assign();

  //// rf pulse should be specified and assigned in seq_block.
  //if (timer_.duration <= 0) {
  //    cout << format("\033[1m\033[31m%s\033[0m %s \033[1m\033[31m'%s'\033[0m %s.\n") % "S-S-L error:" % "RF pulse" % name() % "duration NOT specified yet";
  //    exit(0);
  //}

  //string str;
  //gain_.is_assigned = false;
  //gain_.val = 1;
  //str = getAttribute("beta");
  //if (!str.empty()) {
  //    beta_ = lexical_cast<double>(str);  // unit in deg.
  //    beta_ *= _pi / 180.0;  // into RAD.
  //} else {
  //    // if the flip angle is not specified, we assume the user would like provide
  //    // a 'gain' factor to scale the default pulse amplitude, or give a gain range
  //    // for searching optimal gain factor.
  //    str = getAttribute("gain");
  //    if (!str.empty()) {
  //        vector<double> pars;// = str_cast<double>(str);
  //        if (pars.size() == 1)
  //            gain_.val = pars[0];

  //        // this is specially used when we try to find the optimal gain factor for a fixed
  //        // flip angle beta (e.g. 90 deg) between a given gain interval.
  //        if (pars.size() == 2) {
  //            gain_.min = pars[0];
  //            gain_.max = pars[1];
  //            gain_.is_assigned = true;
  //        }
  //    } else
  //        cout << format("%s %s %s.\n") % "S-S-L warning: " % get_node_name() % "since flip angle has not been specified, RF pulse gain will use default value";
  //}

  //bw_.is_assigned = false;
  //str = getAttribute("bandwidth");
  //if (!str.empty()) {
  //    vector<double> pars;// = str_cast<double>(str);
  //    if (pars.size() == 1)
  //        bw_.val = pars[0];
  //    if (pars.size() > 1) {
  //        bw_.lf = pars[0];
  //        bw_.rf = pars[1];
  //        if (pars.size() > 2)
  //            bw_.pts = pars[2];
  //        else
  //            bw_.pts = 100;
  //        bw_.is_assigned = true;
  //    }
  //}

  //// deal with user defined Tx channels, if not specified in xml node, it will use all Tx
  //// channels in the coil array.
  //str = getAttribute("Tx");
  //if (!str.empty()) {
  //    //tx_index_ = str_cast<size_t>(str);
  //    for (size_t i = 0; i < tx_index_.size(); i++)
  //        tx_index_[i] -= 1;  // switch to 0 based index.
  //}  //else
  ////tx_index_ = pcoil_array_->get_coil_index(kTx);
  //tx_index_.push_back(1);
  //if (tx_index_.size() == 0) {
  //    cout << format("%s %s %s.\n") % "S-S-L error: " % get_node_name() % "block should assign at least 1 Tx coil";
  //    exit(0);
  //}
}

vector<string> rf_pulse::get_channels_str() const {
  vector<string> s;
  for (size_t i = 0; i < raw_data_.size(); i++)
    s.push_back(raw_data_[i].channel);
  return s;
}
int rf_pulse::switch2loop(int index) {
  if (index == -1)
    return 0;
  switch (loop_type_) {
    case ssl::seq::_var_freq: {
      modulated_freq_ = modulated_freq_list_[index];
      cout << modulated_freq_ << "  current freq offset [Hz].\n";
    }
      break;
    case ssl::seq::_var_gain: {
      modulated_gain_ = modulated_gain_list_[index];
      cout << modulated_gain_ << "  current modulated gain.\n";
    }
      break;
    case ssl::seq::_var_phase: {

      if (mode_ == _amp_phase) {
        int m = -32 + index;
        if (m >= 0)
          m = m + 1;

        cout << m << "###\n";

        double delta_k = -_pi / 2;//_pi/2;
        cx_vec shape = raw_data0_[0].envelope;

        for (int i = 0; i < shape.size(); i++)
          shape(i) = cd(shape(i).real(), shape(i).imag() + m * delta_k);

        raw_data_[0].envelope = shape;
      }

    }
      break;
    default:break;
  }
  return 1;
}

void rf_pulse::evolution(int index) {
  if (!g_engine) {
    ssl_color_text("warn", "unknown compute engine, rf evolution ignored.\n");
    return;
  }
  cout << "evolution rf " << index << " " << width_in_ms() << " ms\n";
  switch2loop(index);
  seq_const sc;
  sc.rf_if = 1;
  for (size_t k = 0; k < raw_data_.size(); k++) {
    rf_const tmp;
    tmp.channel = raw_data_[k].channel;
    sc.rf.push_back(tmp);
  }
  for (int i = 1; i < timer_.keys.size(); i++) {
    timeline dt = timer_.keys[i] - timer_.keys[i - 1];
    for (size_t k = 0; k < raw_data_.size(); k++) {
      cd val = raw_data_[k].envelope(i - 1);
      if (mode_ == _amp_phase)
        val = amp2xy(val);
      sc.rf[k].u[cx] = val.real();
      sc.rf[k].u[cy] = val.imag();
      sc.rf[k].u *= modulated_gain_;
      //cout << sc.rf[k].u << "\n";
      sc.rf[k].df = modulated_freq_; // reserved for freq.
    }
    // evolution on current step。
    //auto start = system_clock::now();
    g_engine->evolution(dt, sc);
    //auto end = system_clock::now();
    //auto duration = duration_cast<microseconds>(end - start);
    //cout << "Use Time:" << double(duration.count()) * microseconds::period::num / microseconds::period::den<<" s.\n";
  }
}

void rf_pulse::get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const {
  timeline t0 = key0, t1 = key1;
  if (!adjust_key_interval(t0, t1))
    return;

  if (ctrl.rf.size() != 0)
    ctrl.rf.clear();
  for (size_t k = 0; k < raw_data_.size(); k++) {
    rf_const tmp;
    tmp.channel = raw_data_[k].channel;
    ctrl.rf.push_back(tmp);
  }

  ctrl.rf_if = true;

  int i = (int) (t0 / tl_dt_);
  for (size_t k = 0; k < raw_data_.size(); k++) {
    cd val = raw_data_[k].envelope(i);
    if (mode_ == _amp_phase)
      val = amp2xy(val);
    ctrl.rf[k].u[cx] = val.real();
    ctrl.rf[k].u[cy] = val.imag();
    ctrl.rf[k].u *= modulated_gain_;
    //cout << sc.rf[k].u << "\n";
    //sc.rf[k].df = 0; // reserved for freq.
  }
  //ctrl(amp) = raw_data_[0].envelope(i).real();  // RF magnitude in rad/s.
  //ctrl(phi) = raw_data_[0].envelope(i).imag();  // phase in rad.
  //ctrl(cf) = raw_data_[0].carrier(i);
}

vector<double> rf_pulse::clone_raw_data() const {
  vector<double> vals;
  for (size_t i = 0; i < nsteps_; i++)
    for (size_t j = 0; j < channels_; j++) {
      cd val = raw_data_[j].envelope(i);

      if (mode_ == _ux_uy)
        val *= modulated_gain_;
      else
        val = cd(val.real() * modulated_gain_, val.imag());

      vals.push_back(val.real());
      vals.push_back(val.imag());
    }
  return vals;
}

// \int_0^tp {u^2(t)dt}
double rf_pulse::rf_power() const {
  double p = 0;
  if (mode_ == _amp_phase) {
    for (size_t j = 0; j < channels_; j++)
      p += raw_data_[j].envelope.real().cwiseAbs2().sum();
  }

  if (mode_ == _ux_uy) {
    for (size_t j = 0; j < channels_; j++)
      p += raw_data_[j].envelope.cwiseAbs2().sum();
  }

  p *= dt_ * 1e-6; // Hz^2*s
  //p /= (double) nsteps_;
  //p *= modulated_gain_ * modulated_gain_; // do not forget the amplitude gain (in general 1).
  return p;
}

void rf_pulse::set_gain(double g) {
  modulated_gain_ = g;
}

void rf_pulse::update_raw_data(const double *vals) {
  size_t k = 0;
  for (size_t i = 0; i < nsteps_; i++)
    for (size_t j = 0; j < channels_; j++) {
      raw_data_[j].envelope(i) = cd(vals[k], vals[k + 1]);
      k += 2;
    }
  // gain.
  for (size_t j = 0; j < channels_; j++) {
    if (mode_ == _ux_uy)
      raw_data_[j].envelope /= modulated_gain_;
    else
      raw_data_[j].envelope.real() /= modulated_gain_;
  }
}
vector<cx_vec> rf_pulse::export_signal() const {
  vector<cx_vec> sig;
  for (size_t j = 0; j < channels_; j++)
    sig.push_back(raw_data_[j].envelope);
  return sig;
}
double rf_pulse::sampling_freq() const {
  return 1e3 * double(nsteps_) / width_in_ms();
}

void rf_pulse::switch_rf_mode(string mode) {
  boost::to_lower(mode);
  if (mode == "amp/phase")
    convert2(_amp_phase);
  else if (mode == "ux/uy")
    convert2(_ux_uy);
  else {
    string s = "unknown rf mode: ** " + mode + " ** 'amp/phase' or 'ux/uy' required!";
    throw std::runtime_error(s.c_str());
  }
}

string rf_pulse::get_header() const {
  string s = seq_block::get_header() + "\n";
  s += "# pulse steps: " + to_string(nsteps_) + ".\n";
  if (mode_ == _ux_uy)
    s += "# pulse format: x/y components (Hz).\n";
  else
    s += "# pulse format: amplitude (Hz)/phase (deg).\n";
  s += "# isotope: ";
  for (size_t j = 0; j < channels_; j++) {
    s += raw_data_[j].channel + " ";
  }
  return s;
}

mat rf_pulse::get_shape() const {
  mat m(nsteps_, 2 * channels_);
  for (size_t j = 0; j < channels_; j++) {
    m.col(2 * j) = raw_data_[j].envelope.real();
    m.col(2 * j + 1) = raw_data_[j].envelope.imag();

    if (mode_ == _ux_uy) {
      m.col(2 * j) *= modulated_gain_ / (2 * _pi);
      m.col(2 * j + 1) *= modulated_gain_ / (2 * _pi);
    } else {
      m.col(2 * j) *= modulated_gain_ / (2 * _pi);
      m.col(2 * j + 1) *= 180 / _pi;
    }
  }
  return m;
}

void rf_pulse::h5write(H5File &file, string abbr) const {
  if (abbr.empty())
    abbr = name();
  Group group(file.createGroup(abbr));
  ssl::utility::h5write(file, &group, "header", get_header());
  ssl::utility::h5write(file, &group, "shape", get_shape());
  group.close();
}

void rf_pulse::write(ostream &ostr) const {
  seq_block::write(ostr);
  //ostr << get_header() << "\n";
  ostr << get_shape();
}

void rf_pulse::convert2(envelope_style mode) {
  if (mode_ == mode)
    return;

  if (mode == _ux_uy) {
    for (size_t i = 0; i < raw_data_.size(); i++)
      for (int j = 0; j < raw_data_[i].envelope.size(); j++)
        raw_data_[i].envelope(j) = amp2xy(raw_data_[i].envelope(j));
  }
  if (mode == _amp_phase) {
    for (size_t i = 0; i < raw_data_.size(); i++)
      for (int j = 0; j < raw_data_[i].envelope.size(); j++)
        raw_data_[i].envelope(j) = xy2amp(raw_data_[i].envelope(j)); // rad.
  }
  mode_ = mode;
}
cd amp2xy(cd val) {
  double x = val.real() * cos(val.imag());
  double y = val.real() * sin(val.imag());
  return cd(x, y);
}
cd xy2amp(cd val) {
  return cd(abs(val), deg2rad(phase_in_degree(val)));
}

void rf_pulse::plot() const {
  ofstream ofstr("shape.RF");
  ofstr.precision(3);
  write(ofstr); // already into Hz.
  ofstr.close();

  for (size_t i = 0; i < channels_; i++) {
    string col1 = boost::lexical_cast<string>(2 * i + 1);
    string col2 = boost::lexical_cast<string>(2 * i + 2);
#ifdef GP_SCRIPT_OUTPUT
    ofstream gp("plot_rf.gnu");
#else
    Gnuplot gp;
#endif
    gp << "reset\n";
    gp << terminal_cmd(g_output_terminal);
    string title = name() + " [" + raw_data_[i].channel + "]";
    if (g_output_terminal != "qt")
      gp << "set output "<< "'rf-" << title << "." << g_output_terminal << "'\n";

    gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/xyborder.cfg'\n";
    gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/grid.cfg'\n";
    gp << "load '" << g_project_path << "/share/spin-scenario/config/gnuplot/dark2.pal'\n";

    double dt = width_in_ms() / (double) nsteps_;
    gp << "set xlabel 'Time / ms" << "'\n"; 
    gp << "set style fill solid 0.25\n";
    gp << "set border back\n";
    gp << "set key opaque\n";
    gp << "set key samplen 2\n";

    gp << "set xrange [0:" << width_in_ms() << "]\n";
    if (mode_ == _amp_phase) {
      gp << "set multiplot layout 2,1 title '" << title << "'\n";
      gp << "set ylabel 'Hz'\n";
      //gp << "set ytics 2000\n";
      //gp << "set key top center\n";
      gp << "plot 'shape.RF' u (($0+0.5)*" << dt << "):" << col1 << " title 'amp' with boxes ls 1\n";
      gp << "set ylabel 'deg'\n";
      gp << "set ytics 90\n";
      gp << "set yrange [0:360]\n";
      gp << "plot 'shape.RF' u (($0+0.5)*" << dt << "):" << col2 << " title 'phase' with boxes ls 2\n";
      gp << "unset multiplot\n";
    }
    if (mode_ == _ux_uy) {
      gp << "set title '" << title << "'\n"; // lw 5
      gp << "set ylabel 'amplitude / Hz'\n";
      gp << "plot 'shape.RF' u (($0+0.5)*" << dt << "):" << col1 << " title 'ux' w l ls 1 lw 8, 'shape.RF' u (($0+0.5)*" << dt << "):" << col2 << " title 'uy' w l ls 2 lw 8\n";
    }
#ifdef GP_SCRIPT_OUTPUT
    gp.close();
#endif
  }
}

} /* namespace seq */
} /* namespace ssl */
