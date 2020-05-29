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

#ifndef KERNEL_EQUIP_COIL_FACTORY_H_
#define KERNEL_EQUIP_COIL_FACTORY_H_

#include <kernel/equip/coil.h>
#include <map>

namespace ssl {
namespace equip {

class CoilFactory {
 public:
  CoilFactory();
  virtual ~CoilFactory();
  /// \brief get the coil pointer by key from this factory.
  Coil *get_coil(std::string key) const;

  /// \brief clone a coil object through a xml node.
  ///\warning It will clone a default coil from the factory,
  /// thus only the node name is needed.
  Coil *clone_coil(/*const tinyxml2::XMLElement* node*/) const;

  /// \brief clone a coil object simply by key name.
  ///\warning It will clone a default coil from the factory.
  Coil *clone_coil(std::string key) const;
 private:
  std::map<std::string, Coil *> coils_;
};
} /* namespace equip */
} /* namespace ssl */

#endif /* KERNEL_EQUIP_COIL_FACTORY_H_ */
