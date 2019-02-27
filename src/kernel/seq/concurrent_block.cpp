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

#include "concurrent_block.h"
#include "interval_tree.h"

namespace ssl {
namespace seq {

typedef Interval<seq_block *, timeline> block_interval;
typedef IntervalTree<seq_block *, timeline> block_interval_tree;

concurrent_block::concurrent_block() {
}

concurrent_block::~concurrent_block() {
}

void concurrent_block::assign() {
  if (!is_retrievable("name"))
    config_table_.set("name", "glue block: concurrent");

  style_ = _concurrent;

  glue::assign();
}

void concurrent_block::evolution(int index) {
  allign();
  vector<timeline> keys;
  vector<timeline> terminals;
  for (size_t i = 0; i < sub_blocks_.size(); i++) {

    sub_blocks_[i]->switch2loop(index); // NOTE.
    vector<timeline> sub_keys = sub_blocks_[i]->absolute_keys();
    // Please note grad block only have 2 points: start and end time.
    copy(sub_keys.begin(), sub_keys.end(), back_inserter(keys));

    timeline start = sub_blocks_[i]->start_timeline();
    timeline end = sub_blocks_[i]->end_timeline();
    terminals.push_back(start);
    terminals.push_back(end);
  }

  auto min = std::min_element(std::begin(terminals), std::end(terminals));
  auto max = std::max_element(std::begin(terminals), std::end(terminals));
  timer_.width = *max - *min;

  std::sort(keys.begin(), keys.end());
  auto it = std::unique(keys.begin(), keys.end());
  keys.resize(distance(keys.begin(), it));

  for (size_t k = 0; k < keys.size() - 1; k++) {
    seq_const sc;
    sc.acq.adc = 0;
    sc.rf_if = 0;
    sc.delay_if = 0;
    for (size_t i = 0; i < sub_blocks_.size(); i++)
      sub_blocks_[i]->get_ctrl(keys[k], keys[k + 1], sc);
    g_engine->evolution(keys[k + 1] - keys[k], sc);
  }
  cout << "evolution concurrent-block " << index << " " << width_in_ms() << " ms.\n";
}

void concurrent_block::allign() {
  timeline t0 = sub_blocks_[0]->WIDTH();
  for (size_t i = 1; i < sub_blocks_.size(); i++) {
    string reg = sub_blocks_[i]->align_label();
    // only deal with label allign.
    if (reg != "") {
      timeline ti = sub_blocks_[i]->WIDTH();
      double ms = 0;
      if (reg == "c")
        ms = timeline2ms((t0 - ti)) / 2;
      else if (reg == "l")
        ms = 0;
      else if (reg == "r")
        ms = timeline2ms((t0 - ti));
      else {
        string s = "invalid align parameter: " + reg;
        throw std::runtime_error(s.c_str());
      }
      sub_blocks_[i]->set_align(ms);
    }
  }
  for (size_t i = 0; i < sub_blocks_.size(); i++)
    sub_blocks_[i]->rest_allign(); // in case of cycle concurrent_block.

  if (sub_blocks_[0]->start_timeline() != 0)
    throw std::runtime_error("allign offset of 1st sub-block should be 0.");

  vector<timeline> start_timeline;
  for (size_t i = 1; i < sub_blocks_.size(); i++)
    start_timeline.push_back(sub_blocks_[i]->start_timeline());

  timeline offset = 0;
  auto min = std::min_element(std::begin(start_timeline), std::end(start_timeline));
  if (*min < 0)
    offset = -(*min);
  for (size_t i = 0; i < sub_blocks_.size(); i++)
    sub_blocks_[i]->shift_start_timeline(offset);
}

void concurrent_block::sync() {
  sub_blocks_ = simplify(this);
  vector<int> priors = descending_cycle_priorities();
  if (priors.size() > 1)
    throw std::runtime_error("unique cycle priority required for concurrent blocks.");
  if (priors.size() == 1) {
    set_cycle_priority(priors[0]);
    // determine loop counts.
    for (size_t i = 0; i < sub_blocks_.size(); i++)
      if (sub_blocks_[i]->cycle_priority() != -1) {
        loop_ctrl_.loop_count = sub_blocks_[i]->loop_count();
        break;
      }
  }
}

vector<seq_block *> concurrent_block::simplify(const seq_block *sb) {
  vector<seq_block *> result;
  if (sb->category() == _glue) {
    glue *sb1 = (glue *) sb;
    if (sb1->style() == _concurrent) {
      vector<seq_block *> cur_sub_blocks = sb1->sub_blocks();
      for (size_t i = 0; i < cur_sub_blocks.size(); i++) {
        if (cur_sub_blocks[i]->category() != _glue)
          result.push_back(cur_sub_blocks[i]);
        else {
          vector<seq_block *> tmp = simplify(cur_sub_blocks[i]);
          result.insert(result.end(), tmp.begin(), tmp.end());
        }
      }
    }
  }

  return result;
}

void concurrent_block::write(ostream &ostr) const {
  glue::write(ostr);
}

}
}
