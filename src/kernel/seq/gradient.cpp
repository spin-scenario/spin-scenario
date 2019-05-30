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

#include "gradient.h"
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;

namespace ssl {
namespace seq {

gradient::gradient() {
}

gradient::~gradient() {
}

void gradient::evolution(int index) {
  if (!g_engine) {
    ssl_color_text("warn", "unknown compute engine, delay evolution ignored.\n");
    return;
  }
  switch2loop(index);
  seq_const sc;
  //sc.delay_if = true;
  sc.grad.v[channel_] = grad_area() / width_in_ms();
  g_engine->evolution(WIDTH(), sc);
  cout << "evolution gradient " << index << " " << width_in_ms() << " ms\n";
}

void gradient::assign() {
  config_table_.set("category", _grad);

  string str_channel = retrieve_config_table_str("axis"); // required.
  boost::to_lower(str_channel);
  if (str_channel != "x" && str_channel != "y" && str_channel != "z") {
    string s = "invalid gradient channel: " + str_channel;
    throw std::runtime_error(s.c_str());
  }
  if (str_channel == "x")
    channel_ = _gx;
  if (str_channel == "y")
    channel_ = _gy;
  if (str_channel == "z")
    channel_ = _gz;

  string new_name = retrieve_config_table_str("name");
  config_table_.set("name", new_name + ": " + str_channel);

  // for any gradient, the keys only contain the two boundaries.
  timer_.keys = tlvec(2);
  timer_.keys[0] = 0;
  timer_.keys[1] = timer_.width;

  seq_block::assign();
}

void gradient::plot() const {
  ofstream ofstr("gnu_grad");
  ofstr.precision(4);
  write(ofstr);
  ofstr.close();
  string gp = "plot(\"xlabel<time / ms> ylabel<amplitude / mT/m> gnuplot<unset key> title<" + name() + ">" + "xrange<0:"
      + std::to_string(width_in_ms()) + ">\", line('gnu_grad', 'box'))";
  g_lua->script(gp);
}
}
}
