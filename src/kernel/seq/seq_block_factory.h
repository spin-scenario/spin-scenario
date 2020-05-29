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

#include <kernel/seq/seq_block.h>
#include <memory>
#include <map>

namespace ssl {
namespace seq {

// creates all kinds of seq blocks.
class seq_block_factory {
 public:
  seq_block_factory();
  virtual ~seq_block_factory();
  // get the seq block pointer by key from this factory.
  seq_block *get_seq_block(const std::string key) const;
  // clone a seq block object simply by key name.
  seq_block *clone_seq_block(std::string key) const;
  sol::object clone_seq_block(const sol::table &list, sol::this_state s) const;
 private:
  std::map<std::string, seq_block *> seq_blocks_;
};

} /* namespace seq */
} /* namespace ssl */
