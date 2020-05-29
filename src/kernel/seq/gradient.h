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
#include "seq_block.h"

namespace ssl {
namespace seq {

const std::string grad_func_str[5] = {"phase-encoding", "freq-encoding", "slice-selection", "general-gradient", "gradient list"};
enum grad_func {
  _phase_encoding = 0,
  _freq_encoding,
  _slice_selection,
  _general_grad,
  _list_grads
};

const std::string grad_pattern_str[3] = {"trapezoid", "shaped", "analytical"};
enum grad_pattern {
  _trapezoid = 0,
  _shaped,
  _analytical
};

const char grad_channel_str[3] = {'X', 'Y', 'Z'};
enum grad_channel {
  _gx = 0,
  _gy,
  _gz
};

class gradient : public seq_block {
 public:
  gradient();
  ~gradient();
  virtual gradient *Clone() const = 0;
  virtual void evolution(int index = -1);
  virtual void plot() const;
  virtual double grad_area() const = 0;
 protected:
  virtual void assign();
  grad_channel channel_;
  grad_pattern pattern_;
};

}
}