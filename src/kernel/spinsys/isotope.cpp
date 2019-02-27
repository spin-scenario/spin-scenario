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

#include "isotope.h"
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <kernel/utilities/ssl_config.h>
using namespace boost;

namespace ssl {
namespace spinsys {

vector<nuclear_isotope> isotope::s_isotopes_;

isotope::isotope(const string symbol) : isotope_(NULL) {
  if (!s_isotopes_.size())
    load_nuclear_isotope_database();
  bind_isotope(symbol);
}

isotope::~isotope() {
}

void isotope::bind_isotope(const string symbol) {
  vector<nuclear_isotope>::iterator result = find(s_isotopes_.begin(),
                                                  s_isotopes_.end(),
                                                  symbol);
  if (result == s_isotopes_.end()) {
    string s = "isotope ** " + symbol + " ** is NOT in the current NMR literature!";
    throw std::runtime_error(s.c_str());
  } else
    isotope_ = &(*result);
}
void isotope::load_nuclear_isotope_database() {
  ifstream file;
  string file_name = ssl::utility::g_project_path + "/share/spin-scenario/config/spin.dat";
  file.open(file_name.c_str(), ios::in);
  if (!file) {
    string s = "can't open file " + file_name + " for writing!";
    throw std::runtime_error(s.c_str());
  }
  string line;
  typedef tokenizer<char_separator<char> > tokenizer;
  char_separator<char> sep(" ");  // to be hided.
  while (!file.eof()) {
    getline(file, line);
    tokenizer tok(line, sep);

    if (!starts_with(line, "%")) {
      nuclear_isotope iso;
      int i = 1;
      for (tokenizer::iterator it = tok.begin(); it != tok.end(); i++, ++it) {
        switch (i) {
          case 1:  //protons
            iso.protons = lexical_cast<int>(*it);
            break;
          case 2:  //nucleons
            iso.nucleons = lexical_cast<int>(*it);
            break;
          case 3:  //radioactive *, stable -
            break;
          case 4:  //symbol
            if (iso.nucleons)
              iso.symbol = lexical_cast<string>(iso.nucleons) + (*it);
            else
              iso.symbol = (*it);  // for E
            break;
          case 5:  //name
            iso.name = *it;
            break;
          case 6:  //isotope quantum number
            iso.qn = lexical_cast<float>(*it);
            break;
          case 7:  //gyromagnetic ratio divided by nuclear magneton
            iso.gn = lexical_cast<double>(*it);
            break;
          case 8:  //natural abundance in percent
            break;
          case 9:  //electric quadrupole moment in barn
            break;
          default:break;
        }
      }
      s_isotopes_.push_back(iso);
    }

  }
  file.clear();
  file.close();
}

}
}
