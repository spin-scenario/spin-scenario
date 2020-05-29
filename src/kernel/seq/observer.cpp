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

#include "observer.h"
namespace ssl {
namespace seq {

observer::observer() {
}

observer::~observer() {
}

void observer::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "observer");
  config_table_.set("category", _obser);

  timer_.width = 0;
  config_table_.set("width", timeline2ms(timer_.width)); // unit in ms.
  timer_.keys = tlvec::LinSpaced(1, 0, timer_.width);
  seq_block::assign();
}

void observer::get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const {
}

void observer::evolution(int index) {
  std::cout << "evolution observer " << index << " " << width_in_ms() << " ms\n";
}

}
}