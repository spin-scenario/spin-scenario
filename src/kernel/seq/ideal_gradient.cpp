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

#include "ideal_gradient.h"

namespace ssl {
namespace seq {

ideal_gradient::ideal_gradient() {
  if_amp = false;
  if_area = false;
  if_areas = false;
}

ideal_gradient::~ideal_gradient() {
}

void ideal_gradient::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "ideal gradient"); // the config table has been initialized in seq_block().

  pattern_ = _trapezoid;

  if (is_retrievable("func")) {
    std::string str_func = retrieve_config_table_str("func");
    func_ = _general_grad;
    if (str_func == "phase")
      func_ = _phase_encoding;
    if (str_func == "freq")
      func_ = _freq_encoding;
    if (str_func == "slice")
      func_ = _slice_selection;

    loop_ctrl_.loop_count = 1;
    // trapezoid gradient part.
    switch (func_) {
      case _phase_encoding: {
        double delta_area = 1e3 / (g_seq_param->gamma * g_seq_param->fov[_cy]); // mT/m*ms
        double delta_grad =
            2 * g_seq_param->max_grad / (double) (g_seq_param->matrix[_cy]); // graident increment, mT/m.
        //double delta_grad =
        //2 * g_seq_param->max_grad / (double) (g_seq_param->matrix[_cy] - 1); // graident increment, mT/m.
        double min_duration = delta_area / delta_grad; // ms.
        double duration;
        if (is_retrievable("width")) {
          duration = retrieve_config_table_double("width"); // unit in ms.
          if (duration < min_duration) {
            std::string s = "invalid gradient width for phase encoding, minimum width: "
                + boost::lexical_cast<std::string>(min_duration);
            throw std::runtime_error(s.c_str());
          }
        } else
          config_table_.set("width", min_duration);

        double all_area = delta_area * g_seq_param->matrix[_cy];
        std::string par = boost::lexical_cast<std::string>(-all_area / 2) + ":" + boost::lexical_cast<std::string>(all_area / 2) + ":"
            + boost::lexical_cast<std::string>(g_seq_param->matrix[_cy]);
        //std::string par = "0:" + boost::lexical_cast<std::string>(delta_area*g_seq_param->matrix[_cy]) + ":" + boost::lexical_cast<std::string>(g_seq_param->matrix[_cy]);
        config_table_.set("areas", par);
        loop_ctrl_.loop_count = g_seq_param->matrix[_cy];
      }
        break;
      case _freq_encoding: {
        ////////////////////////NOTE ACQ/////////////////////////////////
        double flat_area = g_seq_param->matrix[_cx] * 1e3 / (g_seq_param->gamma * g_seq_param->fov[_cx]); // mT/m*ms.
        double flat_time = (g_seq_param->matrix[_cx] - 1) / g_seq_param->sw;
        config_table_.set("width", flat_time);
        config_table_.set("area", flat_area);
      }
        break;
      case _slice_selection: {

      }
        break;
      default:break;
    }

  }

  load_shape();

  gradient::assign();
}

void ideal_gradient::get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const {
  timeline t0 = key0, t1 = key1;
  if (!adjust_key_interval(t0, t1))
    return;

  double grad_amp = amp;
  //std::cout<<amp<<"###\n";
  ctrl.grad.v[channel_] = grad_amp;
}

void ideal_gradient::load_shape() {
  // https://github.com/pulseq/pulseq/blob/master/matlab/%2Bmr/makeTrapezoid.m
  // https://github.com/juncy/spin-scenario/blob/v0.1.0/src/kernel/seq/pulse_grad_trapezoid.cpp
  if (is_retrievable("amp")) {
    amp = retrieve_config_table_double("amp"); // unit in mT/m.
    if_amp = true;
    if (fabs(amp) > g_seq_param->max_grad) {
      std::string s = "the given amplitude is beyond limit for this gradient: " + name();
      throw std::runtime_error(s.c_str());
    }
  }
  if (is_retrievable("areas")) {
    std::string par = retrieve_config_table_str("areas"); // unit in mT/m*ms.
    colon_sep val;
    if (parse(par, val)) {

      areas = vec(val.num);
      double delta_s = (val.b - val.a) / (double) val.num;
      for (int i = 0; i < val.num; i++)
        areas[i] = -delta_s * (double) val.num / 2.0 + delta_s * (double) i;
      //areas[i]=delta_s*(double (i)-double(val.num-1)/2);

      //std::cout<<areas<<"\n";
      if_areas = true;
      area = areas[0]; // set 1st phase-encoding area.
    }
  }

  if (is_retrievable("area")) {
    area = retrieve_config_table_double("area"); // unit in mT/m*ms.
    if_area = true;
  }

  if (is_retrievable("width")) {
    double duration = retrieve_config_table_double("width"); // unit in ms.
    timer_.width = ms2timeline(duration);
    duration = timeline2ms(ms2timeline(duration));
    if (!if_amp) {
      amp = area / duration;
    }
  } else {
    std::string s = "must set duration for this gradient: " + name();
    throw std::runtime_error(s.c_str());
  }

  timer_.width = WIDTH();
  //std::cout << timer_.width << "\n";

  tps = tlvec(2);
  tps[0] = 0;
  tps[1] = timer_.width;
}

void ideal_gradient::write(std::ostream &ostr) const {
  seq_block::write(ostr);
  ostr << "# channel: " << grad_channel_str[channel_] << ".\n";
}

void ideal_gradient::plot() const {
}

double ideal_gradient::grad_area() const {
  return if_area ? area : width_in_ms() * amp;
}

int ideal_gradient::switch2loop(int index) {
  int id = seq_block::switch2loop(index);
  if (id == -1)
    return 0;

  if (func_ == _phase_encoding)
    config_table_.set("area", areas[id - 1]);

  load_shape();
  return 1;
}

}
}
