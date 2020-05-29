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
//#include<kernel/utilities/ssl_config.h>
//using namespace ssl::utility;
#include <kernel/physx/engine.h>
using namespace ssl::physx;

namespace ssl {
namespace seq {
/* The sequence blocks are divided into two levels.
Blocks in level 1 include RF block, gradient block, time delay block and signal acquisition block.
Blocks in level 2 inculde two major classs of combination, namely temporal sequential or parallel
assembly of arbitrary combination of all blocks. */
enum block_category {
  _glue = 0,
  _rf,
  _grad,
  _delay,
  _acq,
  _obser
};

struct block_timer {
  timeline width = 0;
  tlvec keys;  // key time points of seq block.
};

enum loop_style {
  _global_cycle,
  _local_array
};
struct loop_ctrl {
  loop_style style = _global_cycle;
  int loop_count = 1;
  int loop_prior = -1; // only for cycle loop.
};

// ctrl for aligned block within concurrent_block.
struct align_ctrl {
  timeline offset_orign = 0;
  timeline offset = 0;
  std::string label = "";
  bool is_specified = false; // flag indicates whether this block has been explicitly aligned via <>.
};

// sequence-block super class.
class seq_block {
 public:
  seq_block();
  virtual ~seq_block();
  virtual seq_block *Clone() const {
    return nullptr;
  }
  // configure this seq-block via lua table.
  void config(const sol::table &t);
  void set_cycle_priority(int val);
  void set_loop_style(loop_style s = _global_cycle);

  void set_align(double ms);
  void set_align(std::string label);
  inline std::string align_label() const {
    return align_ctrl_.label;
  }
  inline bool is_align_specified() const {
    return align_ctrl_.is_specified;
  }

  inline int cycle_priority() const {
    return loop_ctrl_.loop_prior;
  }
  inline int loop_count() const {
    return loop_ctrl_.loop_count;
  }

  inline bool local_loop_if() const {
    if (loop_ctrl_.style == _local_array)
      return true;
    else
      return false;
  }
  /*inline loop_style loop_style() const {
      return loop_ctrl_.style;
  }*/
  virtual void set_name(std::string new_name);
  std::string name() const;
  double width_in_ms() const; // unit in [ms].
  timeline WIDTH() const; // unit in [timeline].
  virtual int switch2loop(int index);
  block_category category() const;
  bool is_retrievable(std::string par) const;
  sol::object retrieve_config_table(std::string key) const;
  std::string retrieve_config_table_str(std::string key) const;
  int retrieve_config_table_int(std::string key) const;
  double retrieve_config_table_double(std::string key) const;
  size_t retrieve_config_table_size_t(std::string key) const;

  inline timeline start_timeline() const {
    return align_ctrl_.offset;
  }
  inline timeline end_timeline() const {
    return align_ctrl_.offset + WIDTH();
  }
  inline void shift_start_timeline(timeline offset) {
    align_ctrl_.offset += offset;
  }
  inline void rest_allign() {
    align_ctrl_.offset = align_ctrl_.offset_orign;
  }

  virtual std::string get_header() const;
  virtual void write(std::ostream &ostr = std::cout) const;
  virtual void plot() const;

  // ONLY for RF.
  virtual void switch_rf_mode(std::string mode) {};
  virtual double rf_power() const {
    return 0;
  }
  //virtual void set_modulated_gain(double g) {};

  // get the control value at the given key time point.
  virtual void get_ctrl(const timeline key0, const timeline key1, seq_const &ctrl) const {};
  virtual std::vector<timeline> absolute_keys() const;
  virtual bool adjust_key_interval(timeline &t0, timeline &t1) const;
  virtual void evolution(int index = -1) {};
  virtual void local_loop_evolution();

  // 944e83cc9b95-4dd9b532f4e9bd0f71c1
  std::string uuid() const;
  // 944e83cc-9b95-4dd9-b532-f4e9bd0f71c1
  std::string uuid_hyphen() const;
  void set_uuid() {
    uuid_ = boost::uuids::random_generator()();
  }
  void create_config_table();
  void copy_config_table();
 protected:
  virtual void assign();
 protected:
  sol::table config_table_; // detail attributes of this seq-block.
  // NOTE: high level cloned blocks will share the unique table with the orginal block.
  // use 'copy_config_table' to re-create a new copy if you wish to deal with the config table further.
  block_timer timer_;
  loop_ctrl loop_ctrl_;
  align_ctrl align_ctrl_;
  boost::uuids::uuid uuid_;
};

class sort_seq_block_cycle_priority {
 public:
  bool operator()(const seq_block *sb1, const seq_block *sb2) {
    return sb1->cycle_priority() > sb2->cycle_priority();
  }
};
} /* namespace seq */
} /* namespace ssl */
