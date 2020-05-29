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

struct phase_acq {
  std::vector<double> deg; // unit in degree.
};

class acquire : public seq_block {
 public:
  acquire();
  ~acquire();
  inline acquire *Clone() const {
    return (new acquire(*this));
  }
  virtual int switch2loop(int index);
  virtual void evolution(int index = -1);
  virtual void get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const;
 protected:
  virtual void assign();
 private:
  size_t npts_;
  double sw_;
  timeline tl_dt_;
  phase_acq loop_phase_list_;
};
} /* namespace seq */
} /* namespace ssl */
