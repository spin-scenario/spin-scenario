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

#include "delay.h"
namespace ssl {
namespace seq {

delay::delay() {

}

delay::~delay() {
}

void delay::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "delay");
  config_table_.set("category", _delay);

  sol::object obj = retrieve_config_table("width"); // required. unit in ms.
  if (obj.get_type() == sol::type::number) {
    loop_delay_list_ = tlvec(1);
    loop_delay_list_[0] = ms2timeline(obj.as<double>());
  }
  if (obj.get_type() == sol::type::string) {
    colon_sep val;
    if (parse(obj.as<std::string>(), val))
      loop_delay_list_ = tlvec::LinSpaced(val.num, ms2timeline(val.a), ms2timeline((val.b)));
  }

  //std::cout << loop_delay_list_ << "\n";

  timer_.width = loop_delay_list_[0];
  timer_.keys = tlvec::LinSpaced(2, 0, timer_.width);

  // determine loop counts.
  loop_ctrl_.loop_count = loop_delay_list_.size();

  seq_block::assign();
}

int delay::switch2loop(int index) {
  int id = seq_block::switch2loop(index);
  if (id == -1)
    return 0;
  timer_.width = loop_delay_list_[id - 1];
  return 1;
}

void delay::evolution(int index) {
  if (!g_engine) {
    ssl_color_text("warn", "unknown compute engine, delay evolution ignored.\n");
    return;
  }
  switch2loop(index);
  seq_const sc;
  sc.delay_if = true;
  g_engine->evolution(WIDTH(), sc);
  std::cout << "evolution delay " << index << " " << width_in_ms() << " ms\n";
}

} /* namespace seq */
} /* namespace ssl */