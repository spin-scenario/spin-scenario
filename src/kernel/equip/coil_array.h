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


#ifndef KERNEL_EQUIP_COIL_ARRAY_H_
#define KERNEL_EQUIP_COIL_ARRAY_H_

#include <kernel/equip/coil_factory.h>

namespace ssl {
namespace equip {

class CoilArray {
 public:
  CoilArray();
  CoilArray(const char *coil_file);
  virtual ~CoilArray();
  /// \brief initialize the input *.xml coil array.
  void load_via_xml(const char *coil_file);
  std::vector<Coil *> get_coil_pointer(CoilMode mode) const;
  std::vector<size_t> get_coil_index(CoilMode mode) const;
  inline std::vector<Coil *> get_coils() const {
    return array_;
  }
  inline size_t get_coils_num() const {
    return array_.size();
  }
  void graphic_view() const;
  std::vector<vec3> get_shape() const;
 private:
  /// \brief generate a coil object corresponding to the xml node.
  void create_coil(/*tinyxml2::XMLElement* node*/);
 private:
  std::vector<Coil *> array_;
  CoilFactory factory_coils_;  ///< all types of coils have been built in this factory.
};

} /* namespace equip */
} /* namespace ssl */

#endif /* KERNEL_EQUIP_COIL_ARRAY_H_ */
