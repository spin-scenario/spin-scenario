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

#include "acquire.h"
#include <chrono>
namespace ssl {
namespace seq {

acquire::acquire() {

}

acquire::~acquire() {
}

void acquire::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "acquisition");
  config_table_.set("category", _acq);
  sol::object npts = retrieve_config_table("np"); // required.
  npts_ = npts.as<size_t>();
  sol::object sw = retrieve_config_table("sw"); // required. unit in Hz.
  sw_ = sw.as<double>();
  sw_ *= 1e-3;
  g_seq_param->sw = sw_; // kHz.

  std::string channels = "1H"; // if not specified, use 1H as default.
  if (is_retrievable("channel")) {
    channels = retrieve_config_table_str("channel"); // required.
    g_seq_param->acq_channel = channels;
  }

  if (is_retrievable("observer")) {
    g_seq_param->observer = retrieve_config_table_str("observer");
//    // may cause a bug when engine is not init.
//    sol::object det = retrieve_config_table("observer");
//    if (g_engine == nullptr)
//      std::cout << "engine not init yet!\n";
//    if (det.get_type() == sol::type::string)
//      g_engine->set_observer(det.as<std::string>());
//    else
//      g_engine->set_observer(det.as<sp_cx_vec>());
  }

  if (is_retrievable("phase")) {
    sol::object phase = retrieve_config_table("phase"); // required.
    std::string phase_par = phase.as<std::string>();

    boost::to_lower(phase_par);

    // ugly code to be rewritten.
    std::vector<char> list;

    for (size_t j = 0; j < phase_par.size(); j++) {
      char s = phase_par.at(j);
      list.push_back(s);
      if (s == 'x' || s == 'y') {
        list.insert(list.end() - 1, ' ');
        list.push_back(' ');
      }
      if (s == '-') {
        list.insert(list.end() - 1, ' ');
        list.push_back(phase_par.at(j + 1));
        list.push_back(' ');
        j += 1;
        continue;
      }
    }
    std::string new_s(&list[0], list.size());
    std::vector<std::string> par_vec;
    boost::trim(new_s);
    boost::split(par_vec, new_s, boost::is_any_of("\t "), boost::token_compress_on);
    // ugly code to be rewritten.

    std::map<std::string, double>::const_iterator iter;
    for (size_t i = 0; i < par_vec.size(); i++) {
      iter = g_phase_map.find(par_vec[i]);
      if (iter != g_phase_map.end())
        loop_phase_list_.deg.push_back(iter->second); // unit in degree.
      else
        loop_phase_list_.deg.push_back(boost::lexical_cast<double>(par_vec[i]));
    }
    if (loop_phase_list_.deg.empty())
      loop_phase_list_.deg.push_back(0); // 'x' phase in default.

    g_seq_param->acq_phase = deg2rad(loop_phase_list_.deg[0]); // update the global acq phase.
    loop_ctrl_.loop_count = loop_phase_list_.deg.size();
  }

  timeline dt = ms2timeline(1.0 / sw_);
  //std::cout << npts_ << "  " << sw_ << "  " << dt << "\n";

  timer_.width = dt * (npts_ - 1);
  config_table_.set("width", timeline2ms(timer_.width)); // unit in ms.
  timer_.keys = tlvec::LinSpaced(npts_, 0, timer_.width);
  tl_dt_ = dt;
  //std::cout << timer_.keys << "\n";
  seq_block::assign();
}

int acquire::switch2loop(int index) {
  if (index == -1)
    return 0;

  int ni = loop_phase_list_.deg.size();
  int idi = (index + 1) % ni;
  if (idi == 0)
    idi = ni;
  double cur_phase = loop_phase_list_.deg[idi - 1];
  g_seq_param->acq_phase = deg2rad(cur_phase); // update the global acq phase.
  std::string s = "phi[" + boost::lexical_cast<std::string>(cur_phase) + "]\n";
  ssl_color_text("seq_phase", s);
  return 1;
}
void acquire::get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const {
  timeline t0 = key0, t1 = key1;
  if (!adjust_key_interval(t0, t1))
    return;

  ctrl.acq.adc = 1;
  ctrl.acq.nps = npts_;
  ctrl.acq.last = false;
  ctrl.acq.index = (int) (t1 / tl_dt_);
  //std::cout << ctrl.acq.index << "----\n";
  if (ctrl.acq.index == (int) (npts_ - 1) || npts_ == 1) {
    ctrl.acq.last = true;
  }
}

void acquire::evolution(int index) {
  if (!g_engine) {
    ssl_color_text("warn", "unknown compute engine, acq evolution ignored.\n");
    return;
  }
  std::cout << "evolution acq " << index << " " << width_in_ms() << " ms\n";
  switch2loop(index);

  if (npts_ == 1) {
    seq_const sc;
    sc.acq.adc = true;
    sc.acq.nps = npts_;
    sc.acq.index = 1;
    sc.acq.last = true;
    g_engine->evolution(0, sc);
  } else
    for (int i = 0; i < timer_.keys.size() - 1; i++) {
      seq_const sc;
      sc.acq.adc = true;
      sc.acq.nps = npts_;
      sc.acq.last = false;
      sc.acq.index = i + 1;
      if (sc.acq.index == (int) (npts_ - 1))
        sc.acq.last = true;

      timeline dt = timer_.keys[i + 1] - timer_.keys[i];
      //auto start = system_clock::now();
      g_engine->evolution(dt, sc);
      //auto end = system_clock::now();
      //auto duration = duration_cast<microseconds>(end - start);
//      std::cout << "Use Time:" << double(duration.count()) * microseconds::period::num / microseconds::period::den
//           << " s.\n";
    }


  // expmv_tspan test.
//  seq_const sc;
//  sc.acq.adc = true;
//  sc.acq.shot = true;
//  sc.acq.nps = npts_;
//  g_engine->evolution(timer_.width, sc);

  //
  //g_engine->reset_rho0();
}
} /* namespace seq */
} /* namespace ssl */
