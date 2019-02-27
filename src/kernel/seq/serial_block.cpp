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

#include "serial_block.h"

namespace ssl {
namespace seq {

serial_block::serial_block() {
}

serial_block::~serial_block() {
}

void serial_block::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "glue block: serial");
  style_ = _serial;
  glue::assign();
}

void serial_block::evolution(int index) {
  if (!g_engine) {
    ssl_color_text("warn", "unknown compute engine, serial_block evolution ignored.\n");
    return;
  }
  vector<int> priors = descending_cycle_priorities();
  imat cycling_matrix = descending_cycle_matrix();
  //cout << cycling_matrix << "\n";
  // #1 #2 #3 ...
  // 0  0  0
  // 1  0  0
  // 2  0  0
  // 0  1  0
  // 1  1  0
  // 2  1  0
  // case 1: all sub-blocks cycle only once.
  if (cycling_matrix.rows() == 0) {
    for (size_t j = 0; j < sub_blocks_.size(); j++) {
      if (sub_blocks_[j]->local_loop_if())
        sub_blocks_[j]->local_loop_evolution();
      else
        sub_blocks_[j]->evolution();
    }
  } else { // case 2: sub-blocks cycles based on the cycling matrix.
    int cycle_num = cycling_matrix.rows();
    for (int i = 0; i < cycle_num; i++) {
      ssl_color_text("info", "global cycling index: " + std::to_string(i + 1) + "/" + std::to_string(cycle_num) + "\n");
      //g_engine->tmp();
      for (size_t j = 0; j < sub_blocks_.size(); j++) {
        int loop_index = -1; // by default, this block is supposed to be non-cycling.
        // retrieve current loop index for this block which denoted as # for cycling block.
        vector<int>::iterator pos = find(priors.begin(), priors.end(), sub_blocks_[j]->cycle_priority());
        if (pos != priors.end()) {
          int prior_id = distance(priors.begin(), pos); // column index of the cycling matrix for this block.
          loop_index = cycling_matrix(i, prior_id);
        }
        // local array loop.
        if (sub_blocks_[j]->local_loop_if())
          sub_blocks_[j]->local_loop_evolution();
        else // global cycling loop.
          sub_blocks_[j]->evolution(loop_index);
      }
      cout << "\n";
    }
  }
}

void serial_block::write(ostream &ostr) const {
  glue::write(ostr);
}
vector<const seq_block *> serial_block::retrieve_cycle_seq_blocks(int prior) const {
  vector<const seq_block *> val;
  for (size_t i = 0; i < sub_blocks_.size(); i++)
    if (sub_blocks_[i]->cycle_priority() == prior)
      val.push_back(sub_blocks_[i]);
  return val;
}

int serial_block::max_cycle_count(int prior) const {
  vector<const seq_block *> sbs = retrieve_cycle_seq_blocks(prior);
  if (sbs.empty())
    return 0;
  vector<int> count;
  for (size_t i = 0; i < sbs.size(); i++)
    count.push_back(sbs[i]->loop_count());
  auto maximum = std::max_element(std::begin(count), std::end(count));
  return *maximum;
}

imat serial_block::descending_cycle_matrix() const {
  vector<int> priors = descending_cycle_priorities();
  vector<int> counts;
  int levels = priors.size();
  for (size_t i = 0; i < priors.size(); i++) {
    int cur_count = max_cycle_count(priors[i]);
    counts.push_back(cur_count);
  }

  vector<vector<int>> cycle_list;
  for (size_t i = 0; i < counts.size(); i++)
    cycle_list = cycle_index_list(cycle_list, counts[i]);

  imat m(cycle_list.size(), levels);
  for (int i = 0; i < m.rows(); i++) {
    Eigen::Map<ivec> v(cycle_list[i].data(), levels);
    //v.reverseInPlace();
    m.row(i) = v.transpose();
  }
  return m;
}

vector<vector<int>> serial_block::cycle_index_list(const vector<vector<int>> &old, int count) const {
  vector<vector<int>> result;
  if (old.size() == 0) {
    for (int i = 0; i < count; i++) {
      vector<int> tmp;
      tmp.push_back(i);
      result.push_back(tmp);
    }
  } else {
    for (int i = 0; i < count; i++) {
      vector<vector<int>> old_copy = old;
      for (size_t j = 0; j < old.size(); j++) {
        old_copy[j].push_back(i);
        result.push_back(old_copy[j]);
      }
    }
  }
  return result;
}
}
}
