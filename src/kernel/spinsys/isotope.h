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

#include <string>
#include <vector>
using namespace std;

namespace ssl {
namespace spinsys {

struct nuclear_isotope {
  int protons;
  int nucleons; // atomic mass in atomic mass unit (protons + neutrons)
  string symbol; // 1H, 13C, ...
  string name; // hydrogen, carbon, ..
  double qn; // spin quantum number
  double gn; // gyromagnetic ratio divided by nuclear magneton
  bool operator==(const string s) {
    return symbol == s;
  }
};
class isotope {
 public:
  isotope(const string symbol = "1H");
  ~isotope();

  inline string symbol() const {
    return isotope_->symbol;
  }
  inline float qn() const {
    return isotope_->qn;
  }
  // get hilbert space size of isotope.
  inline size_t hs() const {
    return 2 * qn() + 1;
  }
  // get gyromagnetic ratio of isotope. 'nuclear magneton' ref to
  // http://en.wikipedia.org/wiki/Gyromagnetic_ratio  & http://en.wikipedia.org/wiki/Nuclear_magneton
  inline double gamma() const {
    double kmu_n = 4.789412272455461e+7; // nuclear magneton
    return isotope_->gn * kmu_n;
  }

  inline bool operator!=(const isotope &I) const {
    return (isotope_ != I.isotope_);
  }
  inline bool operator==(const isotope &I) const {
    return (isotope_ == I.isotope_);
  }
 private:
  void load_nuclear_isotope_database();
  void bind_isotope(const string symbol);
 private:
  nuclear_isotope *isotope_;
  static vector<nuclear_isotope> s_isotopes_;
};

}
}