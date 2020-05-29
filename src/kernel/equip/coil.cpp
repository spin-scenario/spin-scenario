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

#include <kernel/equip/coil.h>

namespace ssl {
namespace equip {

Coil::Coil()
    : id_(0),
      pos_(vec3::Zero()),
      mode_(kTxRx) {
  // TODO Auto-generated constructor stub

}

Coil::~Coil() {
  // TODO Auto-generated destructor stub
}
void Coil::assign() {
  /*std::string str = getAttribute("mode");
  if (!str.empty()) {
      if (str == "Tx")
          mode_ = kTx;
      else if (str == "Rx")
          mode_ = kRx;
      else
          mode_ = kTxRx;
  }*/
}
B1map Coil::create_sensitivity_map_VOI(const ivec3 &dim,
                                       const vec3 &reso,
                                       const vec3 &offset) {
  B1map map;
  //map.B1x =  cube(dim[cy], dim[cx], dim[cz],  fill::zeros);
  //map.B1y =  cube(dim[cy], dim[cx], dim[cz],  fill::zeros);
  // vec3 pos;
  //for (size_t nz = dim[cz] / 2 - 1; nz < dim[cz] / 2 + 1; nz++) {
  //  //for (size_t nz = 0; nz < dim[cz]; nz++) {
  //  pos[cz] = (nz - 0.5 * (dim[cz] - 1)) * reso[cz] + offset[cz];
  //  for (size_t nx = 0; nx < dim[cx]; nx++)
  //    for (size_t ny = 0; ny < dim[cy]; ny++) {
  //      pos[cx] = (nx - 0.5 * (dim[cx] - 1)) * reso[cx] + offset[cx];
  //      pos[cy] = (ny - 0.5 * (dim[cy] - 1)) * reso[cy] + offset[cy];
  //      cd val = get_sensitivity(pos);
  //      map.B1x(ny, nx, nz) = val.real() * cos(val.imag());
  //      map.B1y(ny, nx, nz) = val.real() * sin(val.imag());
  //    }
  //}
  return map;
}
} /* namespace equip */
} /* namespace ssl */
