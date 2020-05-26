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

#include "pulse_seq.h"
#include "rf_pulse.h"
#include "shaped_rf.h"
#include "serial_block.h"
#include "concurrent_block.h"
#include "gradient.h"
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;
#include <chrono>
using namespace chrono;

namespace ssl {
namespace seq {

sol::object multi_shaped_rf(const sol::table &t, sol::this_state s) {
  sol::state_view lua(s);
  sol::table rfs = lua.create_table();

  shaped_rf *rf = new shaped_rf();
  sol::table tt=t;
  tt.set("pattern","rand");
  rf->config(tt);

  string file = retrieve_table_str("file", t);
  hid_t file_id;
  file_id = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  g_h5_string.clear();
  H5Literate(file_id, H5_INDEX_NAME, H5_ITER_NATIVE, NULL, op_func, NULL);
  H5Fclose(file_id);

  H5File h5file;
  h5file.openFile(file, H5F_ACC_RDWR);
  for (size_t i = 0; i < g_h5_string.size(); i++) {
    mat m = h5read_mat(h5file, "/" + g_h5_string[i] + "/shape");
    shaped_rf *p_rf = rf->Clone();
    p_rf->copy_config_table();
    p_rf->set_name(g_h5_string[i]);
    p_rf->set_shape(rf->get_channels_str(), m);
    rfs.add((seq_block *) p_rf);
  }
  h5file.close();
  return rfs;
}

void specgram(const seq_block &rf, const sol::table &t) {
  if (rf.category() != _rf)
    return;
  const rf_pulse &shape = (const rf_pulse &) (rf);
  vector<cx_vec> sig = shape.export_signal();
  // add baseline.
  int n = shape.get_steps();
  double ratio = 0;
  if (is_retrievable("ratio", t))
    ratio = retrieve_table("ratio", t).as<double>();
  int add_n = n * ratio;

  if (ratio>0)
  for (size_t i = 0; i < sig.size(); i++) {
     cx_vec new_sig(n + 2 * add_n);
     new_sig.setZero();
     new_sig.segment(add_n, n) = sig[i];
     sig[i] = new_sig;
  }

  vector<string> channels = shape.get_channels_str();
  double fs = shape.sampling_freq();
  for (size_t i = 0; i < sig.size(); i++) 
    //specgram(sig[i], t, fs, shape.name() + "-" + channels[i]);
    specgram(sig[i], t, fs, channels[i]);
}
void specgram(string file_name, const sol::table &t) {
  mat data = eigen_read(file_name);
  if (data.size() == 0) return;

  cx_vec sig(data.rows());
  sig.setZero();

  string col_str;
  if (is_retrievable("col", t)) {
    col_str = retrieve_table("col", t).as<string>();
    vector<string> cols;
    boost::split(cols, col_str, boost::is_any_of(", "), boost::token_compress_on);
    if (cols.size() == 1)
      sig.real() = data.col(boost::lexical_cast<int>(cols[0]) - 1);
    if (cols.size() >= 2) {
      sig.real() = data.col(boost::lexical_cast<int>(cols[0]) - 1);
      sig.imag() = data.col(boost::lexical_cast<int>(cols[1]) - 1);
    }
  } else
    sig.real() = data.col(0); // by default, use the 1st column of the input file.

  specgram(sig, t, 0);
}
void specgram(const cx_vec &sig, const sol::table &t, double fs, string label) {
  win_shape wshape = utility::_hamming;
  int wlen; // length of the hamming window, recomended to be power of 2.
  int hop; // hop size.
  double overlap;
  int nfft; // number of FFT points, recomended to be power of 2.
  // fs: sampling frequency, Hz

  // deal with the param table.
  if (is_retrievable("window", t))
    wshape = window_interpreter(retrieve_table("window", t).as<string>());

  nfft = retrieve_table("nfft", t, "specgram").as<int>(); // required.
  wlen = retrieve_table("wlen", t, "specgram").as<int>(); // required.
  if (fs == 0)
    fs = retrieve_table("fs", t, "specgram").as<double>(); // required.

  if (is_retrievable("overlap", t))
    overlap = retrieve_table("overlap", t).as<double>();
  else
    overlap = 0.75;

  hop = wlen * (1 - overlap);

  stft_out out = stft(sig, wshape, wlen, hop, nfft, fs);  // sig.normalized()

  string style = "";
  if (is_retrievable("style", t))
    style = retrieve_table("style", t).as<string>();

  int is_ampdB, is_amp, is_phase;
  if (style == "") {  // plot all.
    is_ampdB = 1;
    is_amp = 1;
    is_phase = 1;
  } else {
    vector<string> par_vec;
    boost::split(par_vec, style, boost::is_any_of("\t |"),
                 boost::token_compress_on);
    is_ampdB = std::count(par_vec.begin(), par_vec.end(), "dB");
    is_amp = std::count(par_vec.begin(), par_vec.end(), "amp");
    is_phase = std::count(par_vec.begin(), par_vec.end(), "phase");
  }

  vec2 xrange, yrange;
  xrange[0] = out.time[0];
  xrange[1] = out.time[out.time.size() - 1];
  xrange *= 1e3;
  yrange[0] = out.freq[0];
  yrange[1] = out.freq[out.freq.size() - 1];
  yrange *= 1e-3;
  double f0 = 0;
  double f1 = 0;
  if (is_retrievable("f0", t)) f0 = retrieve_table("f0", t).as<double>();
  if (is_retrievable("f1", t)) f1 = retrieve_table("f1", t).as<double>();

  string frange = "";
  if (f1 > f0)
    frange = "yrange<" + boost::lexical_cast<string>(f0) + ":" + boost::lexical_cast<string>(f1) + ">";

  string fig_spec ="'xlabel<time/ ms> ylabel<freq / kHz> color<Spectral> gnuplot<set palette negative> " +
      frange;
  if (is_amp) {
    utility::map gnu(out.amp);
    gnu.xrange = xrange;
    gnu.yrange = yrange;
    (*g_lua)["_amp"] = gnu;
    g_lua->script("plot(" + fig_spec + " title<magnitude specgram [" + label+ "]>', _amp)");
  }

  if (is_ampdB) {
    utility::map gnu(out.ampdB);
    gnu.xrange = xrange;
    gnu.yrange = yrange;
    (*g_lua)["_ampdB"] = gnu;
    g_lua->script("plot(" + fig_spec +
                  " title<magnitude specgram (dB) [" + label+ "]>', _ampdB)");
  }

  if (is_phase) {
    utility::map gnu(out.unwrap_phase);
    gnu.xrange = xrange;
    gnu.yrange = yrange;
    (*g_lua)["_phase"] = gnu;
    g_lua->script("plot(" + fig_spec + " title<phase specgram [" + label+ "]>', _phase)");
  }
}
void plot(sol::variadic_args va, const seq_block & /*sb*/) {
  for (auto v : va) {
    seq_block &val = v;
    val.plot();
  }
}

void plot(const sol::table &t) {
  for (size_t i = 0; i < t.size(); i++) {
    sol::object val = t[i + 1];
    seq_block *rf = val.as<seq_block *>();
    rf->plot();
  }
}

void write(string file, sol::variadic_args va, const seq_block & /*sb*/) {
  ofstream ofstr(file.c_str());
  ofstr << "# " << sys_time()<<"\n";
  string sep = "#----------------------------------------";
  //ofstr.precision(4);
  for (auto v : va) {
    seq_block &val = v;
    val.write(ofstr);
    ofstr << "\n" << sep;
  }
  ofstr.close();
}
void print(sol::variadic_args va, const seq_block & /*sb*/) {
  string sep = "----------------------------------------\n";
  std::cout << sys_time() << sep;
  std::cout.precision(4);
  for (auto v : va) {
    seq_block &val = v;
    val.write(std::cout);
    std::cout << sep;
  }
}

sol::object run_seq(const seq_block &sb) {
  auto start = system_clock::now();
  cout << "start....\n";

  // in case of non-glue type bolck.
 if (sb.category() != _glue)
    return nullptr;
  // in case of non-serial block.
  glue &serial = (glue &) sb;
  if (serial.style() != _serial)
    return nullptr;
  serial.evolution();

  auto end = system_clock::now();
  auto duration = duration_cast<microseconds>(end - start);
  cout << "Use Time:" << double(duration.count()) * microseconds::period::num / microseconds::period::den << " s.\n";
  cout << "end....\n";

  sol::table result = g_lua->create_table();

  //result.set("raw", double(duration.count()) * microseconds::period::num / microseconds::period::den);
  //return result;//g_engine->process_signal();
  return  g_engine->process_signal();

}

seq_block &serial(const sol::table &t) {
  int num = t.size();
  vector<seq_block *> sbs;
  for (int i = 1; i <= num; i++) {
    sol::object val = t[i];
    switch (val.get_type()) {
      case sol::type::number: {
        std::cout << "number\n";
      }
        break;
      case sol::type::userdata: {
        seq_block &sb = val.as<seq_block &>();
        // deal with the concurrent blocks.
        if (sb.category() == _glue) {
          glue &sb1 = (glue &) sb;
          if (sb1.style() == _concurrent) {
            concurrent_block &sb2 = (concurrent_block &) sb1;
            sb2.sync();
          }
        }
        sbs.push_back(&sb);
        //cout << "running: " << sb.name() << "\n";
        //sb.write(std::cout);
      }
        break;
      case sol::type::table: {
        //std::cout << "table\n";
      }
        break;
      default:break;
    }
  }
  return serial(sbs);
}

seq_block &set_cycle_priority(seq_block &sb, int priority) {
  sb.set_cycle_priority(priority);
  return sb;
}

seq_block &set_loop_array(seq_block &sb, int) {
  sb.set_loop_style(_local_array);
  return sb;
}

seq_block &set_align(seq_block &sb, double ms) {
  sb.set_align(ms);
  return sb;
}
seq_block &set_align(seq_block &sb, string label) {
  sb.set_align(label);
  return sb;
}

seq_block &concurrent(seq_block &sb1, seq_block &sb2) {
  g_lua->script("conc = concurrent{}");
  sol::object val = (*g_lua)["conc"];
  seq_block &sb = val.as<seq_block &>();
  concurrent_block &conc = (concurrent_block &) sb;

  if (!sb2.is_align_specified())
    sb2.set_align("c");

  conc.add_sub_block(&sb1);
  conc.add_sub_block(&sb2);
  return conc;
}

seq_block &serial(vector<seq_block *> sbs) {
  g_lua->script("seri = serial{}");
  sol::object val = (*g_lua)["seri"];
  seq_block &sb = val.as<seq_block &>();
  serial_block &seri = (serial_block &) sb;
  seri.add_sub_block(sbs);
  return seri;
}

engine *init_compute_engine(const sol::table &t) {
  g_engine = new engine(t);
  return g_engine;
}

sol::object run_seq_api(const sol::table &t) {
  init_compute_engine(t);

  sol::object val = retrieve_table("exp", t, "run{}");
  const seq_block &sb = val.as<const seq_block &>();

  sol::table result = g_lua->create_table();

  //result.set("raw", run_seq(sb));
  return run_seq(sb);
}

void init_seq_param(const sol::table &t) {
  if (is_retrievable("fov", t)) {
    sol::object par = retrieve_table("fov", t); // unit in mT/m.
    string par_str = par.as<string>();
    boost::cmatch what;
    boost::regex reg("(\\-?\\d+\\.?\\d*)\\*((\\-?\\d+\\.?\\d*))");
    if (boost::regex_search(par_str.c_str(), what, reg)) {
      double x = boost::lexical_cast<double>(what[1]);
      double y = boost::lexical_cast<double>(what[2]);
      g_seq_param->fov[_cx] = x;
      g_seq_param->fov[_cy] = y;
    }
  }
  if (is_retrievable("matrix", t)) {
    sol::object par = retrieve_table("matrix", t); // unit in mT/m.
    string par_str = par.as<string>();
    boost::cmatch what;
    boost::regex reg("(\\d+)\\*(\\d+)");
    if (boost::regex_search(par_str.c_str(), what, reg)) {
      int x = boost::lexical_cast<int>(what[1]);
      int y = boost::lexical_cast<int>(what[2]);
      g_seq_param->matrix[_cx] = x;
      g_seq_param->matrix[_cy] = y;
    }
  }
  if (is_retrievable("max_grad", t)) {
    sol::object par = retrieve_table("max_grad", t); // unit in mT/m.
    g_seq_param->max_grad = par.as<double>();
  }
  if (is_retrievable("max_slew", t)) {
    sol::object par = retrieve_table("max_slew", t); // unit in T/m/s.
    g_seq_param->max_slew = par.as<double>();
  }

  g_seq_param->write();
}

double area(const seq_block &sb) {
  double s = 0;
  if (sb.category() == _grad) {
    const gradient &sb1 = (const gradient &) sb;
    s = sb1.grad_area();
  } else {
    string s = "area() is only valid for gradient seq-block!";
    throw std::runtime_error(s.c_str());
  }
  return s;
}
}
}
