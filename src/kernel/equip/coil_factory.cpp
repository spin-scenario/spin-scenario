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

#include <kernel/equip/coil_factory.h>
#include <kernel/equip/coil_ideal.h>
#include <kernel/equip/coil_biot_savart.h>
#include <boost/algorithm/string.hpp>

namespace ssl {
namespace equip {

CoilFactory::CoilFactory() {
  // TODO Auto-generated constructor stub
  coils_.insert(std::pair<std::string, Coil *>("COILIDEAL", new CoilIdeal()));
  coils_.insert(std::pair<std::string, Coil *>("COILBIOTSAVART", new CoilBiotSavart()));
}

CoilFactory::~CoilFactory() {
  // TODO Auto-generated destructor stub
}
Coil *CoilFactory::get_coil(std::string key) const {
  boost::to_upper(key);  // in case of case-sensitive key, such as PulseRect, PULSErect, etc.
  return coils_.find(key)->second;
}
Coil *CoilFactory::clone_coil(/*const tinyxml2::XMLElement* node*/) const {
  return nullptr; // clone_coil(/*std_string(node->Name()*/"");
}
Coil *CoilFactory::clone_coil(std::string key) const {
  Coil *to_be_cloned = get_coil(key);
  if (!to_be_cloned)
    return NULL;
  return to_be_cloned->Clone();
}
} /* namespace equip */
} /* namespace ssl */
