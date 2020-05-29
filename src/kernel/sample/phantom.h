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

#pragma once

#include <kernel/utilities/ssl_config.h>
using namespace ssl::utility;

namespace ssl {
namespace sample {

enum phantom_model { mida_brain = 0, mni_brain = 1, usr_phantom =2, unidentified_phantom = 3 };
enum isochromat_prop {
  //cx = 0,
  //cy = 1,
  //cz = 2,
      pd = 3,  // subscript starts at 3, cx, cy, cz already defined in cartesian.
  r1 = 4,
  r2 = 5,
  r2s = 6,
  dB0 = 7
};

// a microscopic group of spins, which resonate at the same frequency.
struct isochromat {
  vec8 data;
  isochromat() {
    data[pd] = 1;
    data[dB0] = 0;
  }
  //std::vector<cd> sens_tx;
  //std::vector<cd> sens_rx;
  // the delta B in x/y/z direction, unit in mT.
  double dB(double gx, double gy, double gz) const {
    double dB = 0;
    if (gx)
      dB += data[cx] * gx;
    if (gy)
      dB += data[cy] * gy;
    if (gz)
      dB += data[cz] * gz;
    return dB;
  }
  // the delta B at this position, unit in mT.
  double delta_B0() const {
    return data[dB0];
  }
  // the spatial location of this voxel, unit in m.
  vec3 position() const {
    vec3 p;
    p(cx) = data[cx];
    p(cy) = data[cy];
    p(cz) = data[cz];
    return p;
  }
};

class phantom {
 public:
  phantom();
  phantom(const char *filename);
  virtual ~phantom();
  void load(const char *filename);
  void init_ensemble();
  void view(const sol::table &t) const;
  inline ivec3 dim() const {
    return dim_;
  }
 private:
  void view(const std::string &axis, int slice) const;
 private:
  ivec3 dim_;  // x,y,z dimension of the cube phantom.
  vec3 res_;   // x,y,z resolution of the cube phantom, unit in m.
  vec3 offset_;

  cube pd_;  // spin density.
  cube T1_;  // T1 time, unit in s.
  cube T2_;  // T2 time, unit in s.
  cube T2s_;
  cube dB0_;

  // ONLY for mida model.
  icube tissue_dist_;
  ivec tissue_index_;
  mat tissue_t1t2_;  //col 2 - T1, col 3 - T2, unit in s.

  phantom_model model_;
 public:
  std::vector<isochromat> isochromats_;
};

}
}
