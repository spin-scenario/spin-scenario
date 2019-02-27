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


#ifndef KERNEL_EQUIP_COIL_BIOT_SAVART_H_
#define KERNEL_EQUIP_COIL_BIOT_SAVART_H_

#include <kernel/equip/coil.h>

namespace ssl {
namespace equip {

class CoilBiotSavart : public Coil {
 public:
  CoilBiotSavart();
  virtual ~CoilBiotSavart();
  inline CoilBiotSavart *Clone() const {
    return (new CoilBiotSavart(*this));
  }
  virtual void assign();
  virtual cd get_sensitivity(const vec3 &pos) const;
 private:
  vector<vec3> get_sampe_points() const;
 private:
  double radius_;  ///< the coil circle radius.
  double azimuth_;  ///< the azimuth angle of the circle plane.
  double elevation_;  ///< the elevation angle of the circle plane.
  size_t segment_;  ///< the number of line segments for approximating circle.
  int current_direction_;  ///< 1 for clockwise, -1 for counterclockwise.
  //vector< vec3> segments_;  ///< segment points on loop.
};

} /* namespace equip */
} /* namespace ssl */

#endif /* KERNEL_EQUIP_COIL_BIOT_SAVART_H_ */
