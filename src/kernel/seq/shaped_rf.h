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
#include "rf_pulse.h"

namespace ssl {
namespace seq {
enum rf_pattern {
  _sinc = 0,
  _gauss,
  _hamming,
  _rand,
  _rect,
  _rand_spline,
  _user_defined
};

std::map<string, rf_pattern> rf_pattern_map();
const std::map<string, rf_pattern> g_rf_pattern = rf_pattern_map();

class shaped_rf : public rf_pulse {
 public:
  shaped_rf();
  virtual ~shaped_rf();
  inline shaped_rf *Clone() const {
    return (new shaped_rf(*this));
  }
  void set_shape(vector<string> channels, rf_pattern p, double maxamp);
  void set_shape(vector<string> channels, string file);
  void set_shape(vector<string> channels, const mat &m);
 protected:
  virtual void assign();

};

} /* namespace seq */
} /* namespace ssl */
