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

#include "analytical_gradient.h"

namespace ssl {
namespace seq {

analytical_gradient::analytical_gradient() {
}

analytical_gradient::~analytical_gradient() {
}

void analytical_gradient::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "expr gradient"); // the config table has been initialized in seq_block().

  pattern_ = _analytical;

  pattern_expr_ = retrieve_config_table_str("expr"); // required.
  if (!pattern_expr_.empty()) {
    pattern_expr_ = yacas_evaluate(pattern_expr_); // in case of D(t) style.
    pattern_integral_header_ = "f" + uuid();
    yacas_evaluate(pattern_integral_header_ + "(t):= (" + yacas_integral(pattern_expr_, "t")
                       + ")"); // need for further key point integral.
  }

  timer_.width = ms2timeline(retrieve_config_table_double("width"));// required.

  gradient::assign();
}

void analytical_gradient::get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const {
  timeline t0 = key0, t1 = key1;
  if (!adjust_key_interval(t0, t1))
    return;

  double grad_amp;
  //std::cout << t0 << "  " << t1 << "\n";
  std::string s =
      "N(" + pattern_integral_header_ + "(" + std::to_string(timeline2ms(t1)) + ")" + "-" + pattern_integral_header_
          + "(" + std::to_string(timeline2ms(t0)) + "))";
  grad_amp = std::stod(yacas_evaluate(s)) / double(timeline2ms((t1 - t0)));
  //std::cout << grad_amp << "\n";

  ctrl.grad.v[channel_] = grad_amp;
}

void analytical_gradient::write(std::ostream &ostr) const {
  seq_block::write(ostr);
  ostr << "# channel: " << grad_channel_str[channel_] << ".\n";
  int np = 256;
  mat raw = mat(np, 2);
  raw.col(0) = vec::LinSpaced(np, 0, width_in_ms());

  yacas_evaluate("f(t):=" + pattern_expr_);
  for (int i = 0; i < np; i++)
    raw(i, 1) = std::stod(yacas_evaluate("N(f(" + std::to_string(raw(i, 0)) + "))"));
  ostr << raw;
}

double analytical_gradient::grad_area() const {
  std::string s =
      "N(" + pattern_integral_header_ + "(" + std::to_string(width_in_ms()) + ") -" + pattern_integral_header_ + "("
          + std::to_string(0) + "))";
  return std::stod(yacas_evaluate(s));
}

int analytical_gradient::switch2loop(int index) {
  return seq_block::switch2loop(index);
}

}
}