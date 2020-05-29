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

#ifndef KERNEL_EQUIP_COIL_H_
#define KERNEL_EQUIP_COIL_H_
#include <kernel/utilities/ssl_config.h>
using namespace ssl::utility;

namespace ssl {
/// \brief equipment package.
namespace equip {

/// \brief Two modes for coil.
enum CoilMode {
  kTx,  ///< RF Transmitting mode.
  kRx,  ///< signal Receiving mode.
  kTxRx  ///< both Tx and Rx mode, the default coil mode.
};

struct B1map {
  //cube B1x;
  //cube B1y;
};

class Coil {
 public:
  Coil();
  virtual ~Coil();
  virtual Coil *Clone() const = 0;
  virtual void assign();
  //virtual double get_sensitivity() const=0;
  virtual cd get_sensitivity(const vec3 &pos) const = 0;
  inline CoilMode get_mode() const {
    return mode_;
  }
  std::vector<vec3> get_shape() const {
    return segments_;
  }
  B1map create_sensitivity_map_VOI(const ivec3 &dim, const vec3 &reso, const vec3 &offset = vec3::Zero());
 protected:
  std::vector<vec3> segments_;  ///< TEMP HERE FOR COIL POS VTK PLOT.
  size_t id_;  ///< coil id.
  vec3 pos_;  ///< center position of this coil.
  CoilMode mode_;  ///< Tx/Rx for RF Transmitting and signal Receiving respectively.
};

} /* namespace equip */
} /* namespace ssl */

#endif /* KERNEL_EQUIP_COIL_H_ */
