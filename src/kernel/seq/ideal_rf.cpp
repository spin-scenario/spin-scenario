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

#include "ideal_rf.h"

namespace ssl {
namespace seq {

ideal_rf::ideal_rf() {
}

ideal_rf::~ideal_rf() {
}

void ideal_rf::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "hard rf pulse"); // the config table has been initialized in seq_block().

  // flip angle.
  sol::object val = retrieve_config_table("beta"); // required. unit in degree.
  beta_ = deg2rad(val.as<double>());
  double ratio = beta_ / (_pi * 0.5);

  double amp90 = _pi * 0.5 / (g_pw90 * 1e-6); // rad.

  nsteps_ = 1;

  double pw; // in ms.
  double amp = amp90;
  if (is_retrievable("width")) {
    pw = retrieve_config_table_double("width"); // required.
    amp = beta_ / (pw * 1e-3); // unit in rad.
  } else
    pw = ratio * g_pw90 * 1e-3;

  timer_.width = ms2timeline(pw); // we suppose the ideal pulse duration to be minimum time (5 us).
  timer_.keys = tlvec::LinSpaced(nsteps_ + 1, 0, timer_.width);
  tl_dt_ = timer_.width;

  string channels = "1H"; // if not specified, use 1H as default.
  if (is_retrievable("channel")) {
    channels = retrieve_config_table_str("channel"); // required.
  }

  vector<string> ch_par_vec;
  boost::split(ch_par_vec, channels, boost::is_any_of("\t, |"), boost::token_compress_on);
  channels_ = ch_par_vec.size();

  // phase.
  if (is_retrievable("phase")) {
    sol::object phase = retrieve_config_table("phase"); // required.
    string phase_par = phase.as<string>();

    boost::to_lower(phase_par);

    vector<string> pha_vec;
    boost::split(pha_vec, phase_par, boost::is_any_of("|"), boost::token_compress_on);

    for (size_t i = 0; i < pha_vec.size(); i++) {
      RFChannel chl;
      chl.envelope = cx_vec::Zero(1);
      chl.carrier = vec::Zero(1);
      chl.channel = ch_par_vec[i];

      phase_channel cur_phase;
      cur_phase.channel = ch_par_vec[i];

      // ugly code to be rewritten.
      vector<char> list;
      for (size_t j = 0; j < pha_vec[i].size(); j++) {
        char s = pha_vec[i].at(j);
        list.push_back(s);
        if (s == 'x' || s == 'y') {
          list.insert(list.end() - 1, ' ');
          list.push_back(' ');
        }
        if (s == '-') {
          list.insert(list.end() - 1, ' ');
          list.push_back(pha_vec[i].at(j + 1));
          list.push_back(' ');
          j += 1;
          continue;
        }
      }
      string new_s(&list[0], list.size());
      vector<string> par_vec;
      boost::trim(new_s);
      boost::split(par_vec, new_s, boost::is_any_of("\t "), boost::token_compress_on);
      // ugly code to be rewritten.

      std::map<string, double>::const_iterator iter;
      for (size_t i = 0; i < par_vec.size(); i++) {
        iter = g_phase_map.find(par_vec[i]);
        if (iter != g_phase_map.end())
          cur_phase.deg.push_back(iter->second); // unit in degree.
        else
          cur_phase.deg.push_back(boost::lexical_cast<double>(par_vec[i]));
      }
      if (cur_phase.deg.empty())
        cur_phase.deg.push_back(0); // 'x' phase in default.

      chl.envelope[0] = cd(amp, cur_phase.deg[0]);
      raw_data_.push_back(chl);
      loop_phase_list_.push_back(cur_phase);
    }

    channels_ = raw_data_.size();
  } else
    for (size_t i = 0; i < channels_; i++) {
      RFChannel chl;
      chl.envelope = cx_vec::Zero(1);
      chl.carrier = vec::Zero(1);
      chl.channel = ch_par_vec[i];

      phase_channel cur_phase;
      cur_phase.channel = ch_par_vec[i];
      cur_phase.deg.push_back(0); // 'x' phase in default.
      loop_phase_list_.push_back(cur_phase);

      chl.envelope[0] = cd(amp, cur_phase.deg[0]);
      raw_data_.push_back(chl);
    }

  for (size_t i = 0; i < raw_data_.size(); i++)
    raw_data_[i].envelope.imag() *= _pi / 180; // phase into rad.

  // determine loop counts.
  vector<int> count;
  for (size_t i = 0; i < loop_phase_list_.size(); i++)
    count.push_back(loop_phase_list_[i].deg.size());
  auto bigest = std::max_element(std::begin(count), std::end(count));
  loop_ctrl_.loop_count = *bigest;

  rf_pulse::assign();
}

int ideal_rf::switch2loop(int index) {
  if (index == -1)
    return 0;

  //cout << index << "\n";

//  int phase0 = 117;
//  //phase0 = deg2rad(phase0);
//
//  if(index==0)
//  tmp_phase.push_back(phase0);
//
//  int cur_phase = tmp_phase.back()+(double)index*phase0;
//  if(index!=0)
//    tmp_phase.push_back(cur_phase);
//
//  cur_phase=cur_phase%360;
//
//  cout<<tmp_phase.back()/phase0<<"###\n";
//
//  cd tmp = raw_data_[0].envelope[0];
//  if (mode_ == _ux_uy)
//    tmp = xy2amp(raw_data_[0].envelope[0]);
//  tmp = cd(tmp.real(), deg2rad(double(cur_phase)));
//
//  raw_data_[0].envelope[0] = tmp;
//  if (mode_ == _ux_uy)
//    raw_data_[0].envelope[0] = amp2xy(tmp);


  for (size_t i = 0; i < loop_phase_list_.size(); i++) {
    int ni = loop_phase_list_[i].deg.size();
    int idi = (index + 1) % ni;
    if (idi == 0)
      idi = ni;
    double cur_phase = loop_phase_list_[i].deg[idi - 1];
    string s = "phi[" + boost::lexical_cast<string>(cur_phase) + "]\n";
    ssl_color_text("seq_phase", s);
    cur_phase = deg2rad(cur_phase);
    cd tmp = raw_data_[i].envelope[0];
    if (mode_ == _ux_uy)
      tmp = xy2amp(raw_data_[i].envelope[0]);
    tmp = cd(tmp.real(), cur_phase);

    raw_data_[i].envelope[0] = tmp;
    if (mode_ == _ux_uy)
      raw_data_[i].envelope[0] = amp2xy(tmp);
  }
  return 1;
}

}
}
