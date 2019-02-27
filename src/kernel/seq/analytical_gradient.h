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
#include "gradient.h"

namespace ssl {
namespace seq {

class analytical_gradient : public gradient {
 public:
  analytical_gradient();
  ~analytical_gradient();
  inline analytical_gradient *Clone() const {
    return (new analytical_gradient(*this));
  }
  virtual void get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const;
  virtual void write(ostream &ostr = cout) const;
  virtual double grad_area() const;
  virtual int switch2loop(int index);
 protected:
  virtual void assign();
 private:
  string pattern_expr_;
  string pattern_integral_header_;
};

}
}