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

#include "trapezoid_gradient.h"

namespace ssl {
namespace seq {

trapezoid_gradient::trapezoid_gradient() {
  if_rise_time = false;
  if_amp = false;
  if_flat_area = false;
  if_area = false;
  if_areas = false;
}

trapezoid_gradient::~trapezoid_gradient() {
}

void trapezoid_gradient::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "trapezoid gradient"); // the config table has been initialized in seq_block().

  pattern_ = _trapezoid;

  if (is_retrievable("func")) {
    string str_func = retrieve_config_table_str("func");
    func_ = _general_grad;
    if (str_func == "phase_encode")
      func_ = _phase_encoding;
    if (str_func == "read_out")
      func_ = _freq_encoding;
    if (str_func == "slice_select")
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
            string s = "invalid gradient width for phase encoding, minimum width: "
                + boost::lexical_cast<string>(min_duration);
            throw std::runtime_error(s.c_str());
          }
        } else
          config_table_.set("width", min_duration);

        double all_area = delta_area * g_seq_param->matrix[_cy];
        string par = boost::lexical_cast<string>(-all_area / 2) + ":" + boost::lexical_cast<string>(all_area / 2) + ":"
            + boost::lexical_cast<string>(g_seq_param->matrix[_cy]);
        //string par = "0:" + boost::lexical_cast<string>(delta_area*g_seq_param->matrix[_cy]) + ":" + boost::lexical_cast<string>(g_seq_param->matrix[_cy]);
        config_table_.set("areas", par);
        loop_ctrl_.loop_count = g_seq_param->matrix[_cy];
      }
        break;
      case _freq_encoding: {
        ////////////////////////NOTE ACQ/////////////////////////////////
        double flat_area = g_seq_param->matrix[_cx] * 1e3 / (g_seq_param->gamma * g_seq_param->fov[_cx]); // mT/m*ms.
        double flat_time = (g_seq_param->matrix[_cx] - 1) / g_seq_param->sw;
        config_table_.set("flat_time", flat_time);
        config_table_.set("flat_area", flat_area);
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

void trapezoid_gradient::get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const {
  timeline t0 = key0, t1 = key1;
  if (!adjust_key_interval(t0, t1))
    return;

  double grad_amp = 0;
  vector<grad_interval> slected_sub_intervals;
  tree.findOverlapping(t0, t1, slected_sub_intervals);
  if (slected_sub_intervals.size() == 1) {
    string s = "N(" + slected_sub_intervals[0].value.header + "(" + std::to_string((t1)) + ") -"
        + slected_sub_intervals[0].value.header + "(" + std::to_string((t0)) + "))";
    grad_amp = std::stod(yacas_evaluate(s)) / double(t1 - t0);
  } else if (slected_sub_intervals.size() == 2) {
    timeline tk = tps[slected_sub_intervals[1].value.pos];
    string s1 = "N(" + slected_sub_intervals[0].value.header + "(" + std::to_string((tk)) + ") -"
        + slected_sub_intervals[0].value.header + "(" + std::to_string((t0)) + "))";
    string s2 = "N(" + slected_sub_intervals[1].value.header + "(" + std::to_string((t1)) + ") -"
        + slected_sub_intervals[1].value.header + "(" + std::to_string((tk)) + "))";
    grad_amp = (std::stod(yacas_evaluate(s1)) + std::stod(yacas_evaluate(s2))) / double(t1 - t0);
  } else if (slected_sub_intervals.size() == 3) {
    string s1 = "N(" + slected_sub_intervals[0].value.header + "(" + std::to_string((tps[1])) + ") -"
        + slected_sub_intervals[0].value.header + "(" + std::to_string((t0)) + "))";
    string s2 = "N(" + slected_sub_intervals[1].value.header + "(" + std::to_string((tps[2])) + ") -"
        + slected_sub_intervals[1].value.header + "(" + std::to_string((tps[1])) + "))";
    string s3 = "N(" + slected_sub_intervals[2].value.header + "(" + std::to_string((t1)) + ") -"
        + slected_sub_intervals[2].value.header + "(" + std::to_string((tps[2])) + "))";
    grad_amp = (std::stod(yacas_evaluate(s1)) + std::stod(yacas_evaluate(s2)) + std::stod(yacas_evaluate(s3)))
        / double(t1 - t0);
  }
  //cout << "[" << t0 << " " << t1 << "]" << "  " << grad_amp << "\n";
  ctrl.grad.v[channel_] = grad_amp;
}

void trapezoid_gradient::load_shape() {
  // https://github.com/pulseq/pulseq/blob/master/matlab/%2Bmr/makeTrapezoid.m
  // https://github.com/juncy/spin-scenario/blob/v0.1.0/src/kernel/seq/pulse_grad_trapezoid.cpp
  if (is_retrievable("amp")) {
    amp = retrieve_config_table_double("amp"); // unit in mT/m.
    if_amp = true;
    if (fabs(amp) > g_seq_param->max_grad) {
      string s = "the given amplitude is beyond limit for this gradient: " + name();
      throw std::runtime_error(s.c_str());
    }
  }

  if (is_retrievable("flat_area")) {
    flat_area = retrieve_config_table_double("flat_area"); // unit in mT/m*ms.
    if_flat_area = true;
  }

  if (is_retrievable("areas")) {
    string par = retrieve_config_table_str("areas"); // unit in mT/m*ms.
    colon_sep val;
    if (parse(par, val)) {

      areas = vec(val.num);
      double delta_s = (val.b - val.a) / (double) val.num;
      for (int i = 0; i < val.num; i++)
        //areas[i]=-delta_s*(double)val.num/2.0+delta_s*(double)i;
        areas[i] = delta_s * (i - val.num / 2);

      //cout<<areas<<"\n";
      if_areas = true;
      area = areas[0]; // set 1st phase-encoding area.
    }
  }

  if (is_retrievable("area")) {
    area = retrieve_config_table_double("area"); // unit in mT/m*ms.
    if_area = true;
  }

  if (is_retrievable("rise_time")) {
    rise_time = retrieve_config_table_double("rise_time"); // unit in ms.
    if_rise_time = true;
  }

  // flat_time + amp or flat_area
  // duration + amp or area
  if (is_retrievable("flat_time")) {
    flat_time = retrieve_config_table_double("flat_time"); // unit in ms.
    flat_time = timeline2ms(ms2timeline(flat_time));
    if (!if_amp)
      amp = flat_area / flat_time;

    if (!if_rise_time)
      rise_time = fabs(amp) / g_seq_param->max_slew;

    rise_time = timeline2ms(ms2timeline(rise_time));
  } else if (is_retrievable("width")) {
    double duration = retrieve_config_table_double("width"); // unit in ms.
    duration = timeline2ms(ms2timeline(duration));
    if (!if_amp) {
      bool possible;
      if (if_rise_time) {
        amp = area / (duration - rise_time);
        possible = fabs(amp) <= g_seq_param->max_grad && (duration - 2 * rise_time) > 0;
      } else {
        double dC = 1 / abs(2 * g_seq_param->max_slew) + 1 / abs(2 * g_seq_param->max_slew);
        amp = (duration - sqrt(duration * duration - 4 * abs(area) * dC)) / (2 * dC);
        if (area < 0)
          amp *= -1;
        possible = fabs(amp) <= g_seq_param->max_grad && (duration * duration > 4 * abs(area) * dC);
      }
      if (!possible) {
        string s = "the given area is too large for this gradient: " + name();
        throw std::runtime_error(s.c_str());
      }
    }
    if (!if_rise_time)
      rise_time = fabs(amp) / g_seq_param->max_slew;

    rise_time = timeline2ms(ms2timeline(rise_time));
    flat_time = timeline2ms(ms2timeline((duration - 2 * rise_time)));
    //cout << rise_time<<" " << flat_time << "\n";
    // adjust amplitude (after rounding) to achieve given area.
    if(!if_amp)
    amp = area / (flat_time + rise_time);
  } else {
    string s = "must set duration for this gradient: " + name();
    throw std::runtime_error(s.c_str());
  }

  timer_.width = ms2timeline((rise_time * 2 + flat_time));
  //cout << timer_.width << "\n";

  tps = tlvec(4);
  tps[0] = 0;
  tps[1] = ms2timeline(rise_time);
  tps[2] = ms2timeline((rise_time + flat_time));
  tps[3] = ms2timeline((rise_time * 2 + flat_time));

  pattern_expr[0] = std::to_string(amp / tps[1]) + "*x";
  pattern_expr[1] = std::to_string(amp);
  pattern_expr[2] = std::to_string(-amp / tps[1]) + "*(x-" + std::to_string(tps[3]) + ")";

  pattern_integral_header[0] = "f1" + uuid();
  pattern_integral_header[1] = "f2" + uuid();
  pattern_integral_header[2] = "f3" + uuid();
  for (int i = 0; i < 3; i++) {
    pattern_integral_expr[i] = yacas_integral(pattern_expr[i]);
    yacas_evaluate(pattern_integral_header[i] + "(x):= (" + pattern_integral_expr[i] + ")");
  }

  sub_intervals.clear(); // Note this function will be done multi-time, so each time it must be cleared.
  sub_intervals.push_back(grad_interval(tps[0], tps[1], aux_grad_tree(pattern_integral_header[0], 0)));
  sub_intervals.push_back(grad_interval(tps[1], tps[2], aux_grad_tree(pattern_integral_header[1], 1)));
  sub_intervals.push_back(grad_interval(tps[2], tps[3], aux_grad_tree(pattern_integral_header[2], 2)));
  tree = grad_interval_tree(sub_intervals);
}

void trapezoid_gradient::write(ostream &ostr) const {
  seq_block::write(ostr);
  ostr << "# channel: " << grad_channel_str[channel_] << ".\n";
  ostr << 0 << " " << 0 << "\n";
  ostr << rise_time << " " << amp << "\n";
  ostr << rise_time + flat_time << " " << amp << "\n";
  ostr << rise_time * 2 + flat_time << " " << 0 << "\n";
}

trap_data trapezoid_gradient::switch2loop(int index) const {
  double new_area = areas[index];
  trap_data data;
  if (is_retrievable("width")) {
    double duration = retrieve_config_table_double("width"); // unit in ms.
    duration = timeline2ms(ms2timeline(duration));
    if (!if_amp) {
      bool possible;
      if (if_rise_time) {
        data.amp = new_area / (duration - rise_time);
        possible = fabs(amp) <= g_seq_param->max_grad && (duration - 2 * rise_time) > 0;
      } else {
        double dC = 1 / abs(2 * g_seq_param->max_slew) + 1 / abs(2 * g_seq_param->max_slew);
        data.amp = (duration - sqrt(duration * duration - 4 * abs(new_area) * dC)) / (2 * dC);
        if (new_area < 0)
          data.amp *= -1;
        possible = fabs(amp) <= g_seq_param->max_grad && (duration * duration > 4 * abs(new_area) * dC);
      }
      if (!possible) {
        string s = "the given area is too large for this gradient: " + name();
        throw std::runtime_error(s.c_str());
      }
    }
    if (!if_rise_time)
      data.rise_time = fabs(amp) / g_seq_param->max_slew;

    data.rise_time = timeline2ms(ms2timeline(rise_time));
    data.flat_time = timeline2ms(ms2timeline((duration - 2 * rise_time)));
    // adjust amplitude (after rounding) to achieve given area.
    data.amp = new_area / (flat_time + rise_time);
  }

  return data;
}

void trapezoid_gradient::plot() const {
  if (func_ == _phase_encoding) {
    g_lua->script("os.execute('rm -rf grad')");
    g_lua->script("os.execute('mkdir grad')");
    for (int i = 0; i < areas.size(); i++) {
      trap_data data = switch2loop(i);
      ofstream ofstr("grad/phase-" + std::to_string(i + 1));
      ofstr.precision(4);
      ofstr << 0 << " " << 0 << "\n";
      ofstr << data.rise_time << " " << data.amp << "\n";
      ofstr << data.rise_time + data.flat_time << " " << data.amp << "\n";
      ofstr << data.rise_time * 2 + data.flat_time << " " << 0 << "\n";
      ofstr.close();
    }
    string gp =
        "plot(\"xlabel<time / ms> ylabel<amplitude / mT/s> gnuplot<unset key> title<" + grad_func_str[func_] + ">"
            + "xrange<0:" + std::to_string(width_in_ms()) + ">\", lines('grad/phase-1:"
            + std::to_string(areas.size()) + "'))";
    // gnuplot[set key outside]
    g_lua->script(gp);
    return;
  }

  gradient::plot();
}

double trapezoid_gradient::grad_area() const {
  return if_area ? area : (rise_time + flat_time) * amp;
}

int trapezoid_gradient::switch2loop(int index) {
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
