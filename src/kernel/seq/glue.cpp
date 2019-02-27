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

#include "glue.h"
#include <algorithm>

namespace ssl {
namespace seq {
glue::glue() {
}

glue::~glue() {
  // TODO Auto-generated destructor stub
}

void glue::assign() {
  config_table_.set("category", _glue);
  seq_block::assign();
}

void glue::add_sub_block(seq_block *sb) {
  sub_blocks_.push_back(sb);
}
void glue::add_sub_block(vector<seq_block *> sbs) {
  //for each (auto var in sbs) {
  //    sub_blocks_.push_back(var);
  //}
  for (size_t i = 0; i < sbs.size(); i++)
    sub_blocks_.push_back(sbs[i]);
}
void glue::add_sub_block(sol::variadic_args va, const seq_block & /*sb*/) {
  for (auto v : va) {
    seq_block &val = v;
    sub_blocks_.push_back(&val);
  }
}

void glue::write(ostream &ostr) const {
  seq_block::write(ostr);
  for (size_t i = 0; i < sub_blocks_.size(); i++)
    ostr << sub_blocks_[i]->name() << "\n";
}

vector<int> glue::descending_cycle_priorities() const {
  vector<seq_block *> sub_blocks = sub_blocks_;
  std::sort(sub_blocks.begin(), sub_blocks.end(), sort_seq_block_cycle_priority());
  vector<int> val;
  for (size_t i = 1; i < sub_blocks.size(); i++) {
    int cur = sub_blocks[i]->cycle_priority();
    int pre = sub_blocks[i - 1]->cycle_priority();
    if (cur != pre)
      val.push_back(pre);
  }
  int last = sub_blocks.back()->cycle_priority();
  if (last != -1)
    val.push_back(last);

  std::reverse(val.begin(), val.end());
  return val;
}

}
}
