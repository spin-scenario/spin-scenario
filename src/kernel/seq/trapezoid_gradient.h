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
#include "interval_tree.h"

namespace ssl {
namespace seq {

struct aux_grad_tree {
  std::string header;
  int pos;
  aux_grad_tree(std::string a, int b) {
    header = a;
    pos = b;
  }
};
typedef Interval<aux_grad_tree, timeline> grad_interval;
typedef IntervalTree<aux_grad_tree, timeline> grad_interval_tree;

struct trap_data {
  double rise_time; // ms.
  double flat_time; // ms.
  double amp; // mT/m.
};

class trapezoid_gradient : public gradient {
 public:
  trapezoid_gradient();
  ~trapezoid_gradient();
  inline trapezoid_gradient *Clone() const {
    return (new trapezoid_gradient(*this));
  }
  virtual void get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const;
  virtual void write(std::ostream &ostr = std::cout) const;
  virtual void plot() const;
  virtual double grad_area() const;
  virtual int switch2loop(int index);
 protected:
  virtual void assign();
  void load_shape();
  trap_data switch2loop(int index) const;
 private:
  double rise_time; // ms.
  double flat_time; // ms.
  double amp; // mT/m.
  double flat_area; // mT/m*ms.
  double area; // mT/m*ms.
  vec areas;
  bool if_rise_time;
  bool if_amp;
  bool if_flat_area;
  bool if_area;
  bool if_areas;
  grad_func func_;

  tlvec tps;
  std::string pattern_expr[3];
  std::string pattern_integral_expr[3];
  std::string pattern_integral_header[3];
  std::vector<grad_interval> sub_intervals;
  grad_interval_tree tree;
};

}
}
