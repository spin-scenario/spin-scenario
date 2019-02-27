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

#include "seq_block_factory.h"
#include "shaped_rf.h"
#include "ideal_rf.h"
#include "delay.h"
#include "acquire.h"
#include "concurrent_block.h"
#include "serial_block.h"
#include "analytical_gradient.h"
#include "trapezoid_gradient.h"
#include "ideal_gradient.h"
#include "observer.h"

namespace ssl {
namespace seq {

seq_block_factory::seq_block_factory() {
  // low level (L1) blocks, including time delay, RF pulses,
  // gradient pulses and signal acquisition unit.
  seq_blocks_.insert(pair<string, seq_block *>("delay", new delay()));

  /*seq_blocks_.insert(pair<string, seq_block*>("RFRECT", new PulseRFRect()));
  seq_blocks_.insert(pair<string, seq_block*>("RFSINC", new PulseRFSinc()));
  seq_blocks_.insert(pair<string, seq_block*>("RFSLR", new PulseRFSLR()));
  seq_blocks_.insert(pair<string, seq_block*>("RFGAUSSIAN", new PulseRFGaussian()));*/
  seq_blocks_.insert(pair<string, seq_block *>("shapedrf", new shaped_rf()));
  seq_blocks_.insert(pair<string, seq_block *>("idealrf", new ideal_rf()));

  seq_blocks_.insert(pair<string, seq_block *>("analyticalgrad", new analytical_gradient()));
  seq_blocks_.insert(pair<string, seq_block *>("trapezoidgrad", new trapezoid_gradient()));
  seq_blocks_.insert(pair<string, seq_block *>("idealgrad", new ideal_gradient()));
  /*seq_blocks_.insert(pair<string, seq_block*>("GRADTRAP", new PulseGradTrapezoid()));
  seq_blocks_.insert(pair<string, seq_block*>("GRADTRIA", new PulseGradTriangle()));
  seq_blocks_.insert(pair<string, seq_block*>("GRADSPIR", new PulseGradSpiral()));*/

  seq_blocks_.insert(pair<string, seq_block *>("acquire", new acquire()));

  // high level (L2) blocks.
  seq_blocks_.insert(pair<string, seq_block *>("concurrent", new concurrent_block()));
  seq_blocks_.insert(pair<string, seq_block *>("serial", new serial_block()));
  //seq_blocks_.insert(pair<string, seq_block*>("SEQCONFIG", new SeqConfig()));

  seq_blocks_.insert(pair<string, seq_block *>("observer", new observer()));
}

seq_block_factory::~seq_block_factory() {
}

seq_block *seq_block_factory::get_seq_block(const string key) const {
  string upper_key = boost::to_lower_copy(key);  // in case of case-sensitive key, such as PulseRect, PULSErect, etc.
  map<string, seq_block *>::const_iterator iter;
  iter = seq_blocks_.find(upper_key);
  if (iter != seq_blocks_.end())
    return iter->second;
  else {
    string s = "invalid seq_block type ** " + key + " **";
    throw std::runtime_error(s.c_str());
  }
}

seq_block *seq_block_factory::clone_seq_block(string key) const {
  seq_block *to_be_cloned = get_seq_block(key);
  return to_be_cloned->Clone();
}

sol::object seq_block_factory::clone_seq_block(const sol::table &list, sol::this_state s) const {
  sol::state_view lua(s);
  sol::table t = lua.create_table();
  vector<string> keys;
  string last_key;
  for (auto &kv : list) {
    sol::object val = kv.second;
    switch (val.get_type()) {
      case sol::type::number: {
        int num = val.as<int>();
        if (last_key.empty()) {
#ifdef SSL_OUTPUT_ENABLE
          string s =
              str(boost::format("%s %s (%s).\n") % "in seq_block{...}, no seq block name defined before the number:"
                      % num % "ignored");
          ssl_color_text("warn", s);
#endif
          break;
        }
        if (num <= 0) {
#ifdef SSL_OUTPUT_ENABLE
          string s = str(boost::format("%s %s (%s %s).\n")
                             % "in seq_block{...}, each seq block size should be >=1, invalid number:" % num % last_key
                             % "ignored");
          ssl_color_text("warn", s);
#endif
          keys.pop_back();
        }
        for (int i = 1; i < num; i++)
          keys.push_back(last_key);
      }
        break;
      case sol::type::string: {
        string key = val.as<string>();
        keys.push_back(key);
        last_key = key;
      }
        break;
      case sol::type::table:break;
      default:break;
    }
  }

  for (size_t i = 0; i < keys.size(); i++) {
    seq_block *to_be_cloned = get_seq_block(keys[i]);
    t.add((to_be_cloned->Clone()));
  }
  return t;
}
} /* namespace seq */
} /* namespace ssl */
