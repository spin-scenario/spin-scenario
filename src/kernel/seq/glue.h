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

enum glue_style {
  _concurrent,
  _serial
};
class glue : public seq_block {
 public:
  glue();
  virtual ~glue();
  virtual glue *Clone() const = 0;
  virtual void add_sub_block(seq_block *sb);
  virtual void add_sub_block(std::vector<seq_block *> sbs);
  virtual void add_sub_block(sol::variadic_args va, const seq_block & /*sb*/);
  virtual void write(std::ostream &ostr = std::cout) const;
  inline std::vector<seq_block *> sub_blocks() const {
    return sub_blocks_;
  }
  inline glue_style style() const {
    return style_;
  }
 protected:
  virtual void assign();
  // all cycle priorities in descending order (from high to low).
  std::vector<int> descending_cycle_priorities() const;
 protected:
  std::vector<seq_block *> sub_blocks_;
  glue_style style_;
};
}
}