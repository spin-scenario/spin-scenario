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

struct phase_channel {
  std::string channel;
  std::vector<double> deg; // unit in degree.
};

class ideal_rf : public rf_pulse {
 public:
  ideal_rf();
  ~ideal_rf();
  inline ideal_rf *Clone() const {
    return (new ideal_rf(*this));
  }
  virtual int switch2loop(int index);
 protected:
  virtual void assign();
 protected:
  double beta_;  ///< flip angle, unit in rad.
  std::vector<phase_channel> loop_phase_list_;
  std::vector<int> tmp_phase;
};

}
}