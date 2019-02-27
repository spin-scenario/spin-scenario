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

#include "seq_block.h"

namespace ssl {
namespace seq {

seq_param *g_seq_param = new seq_param;

seq_block::seq_block() {
}

seq_block::~seq_block() {
  // TODO Auto-generated destructor stub
}

void seq_block::create_config_table() {
  if (!config_table_.valid())
    config_table_ = g_lua->create_table();
}

void seq_block::copy_config_table() {
  sol::table new_table = g_lua->create_table();
  for (auto &kv : config_table_) {
    sol::object key = kv.first;
    sol::object val = kv.second;
    new_table.set(key, val); // existed keys will be overwrited.
  }
  config_table_ = new_table;
}

string seq_block::uuid() const {
  string s = uuid_hyphen();
  boost::erase_all(s, "-");
  return s;
}
string seq_block::uuid_hyphen() const {
  return boost::lexical_cast<string>(uuid_);
}

void seq_block::config(const sol::table &t) {
  /*if (!t.valid() || t.empty())
      throw std::exception("invalid seq_block parameter table (nil or empty).");*/
  create_config_table();
  for (auto &kv : t) {
    sol::object key = kv.first;
    sol::object val = kv.second;
    config_table_.set(key, val); // existed keys will be overwrited.
  }

  // in case of duplicate uuids of each block type for clone reason.
  uuid_ = boost::uuids::random_generator()();

  assign();
}
void seq_block::set_cycle_priority(int val) {
  loop_ctrl_.loop_prior = val;
}
void seq_block::set_loop_style(loop_style s) {
  loop_ctrl_.style = s;
}

void seq_block::set_align(double ms) {
  align_ctrl_.offset_orign = ms2timeline(ms);
  align_ctrl_.is_specified = true;
}
void seq_block::set_align(string label) {
  align_ctrl_.label = label;
  align_ctrl_.is_specified = true;
}

bool seq_block::is_retrievable(string par) const {
  return ssl::utility::is_retrievable(par, config_table_);
}

sol::object seq_block::retrieve_config_table(string key) const {
  return retrieve_table(key, config_table_, "seq_block config");
}

string seq_block::retrieve_config_table_str(string key) const {
  return retrieve_table_str(key, config_table_, "seq_block config");
}

int seq_block::retrieve_config_table_int(string key) const {
  return retrieve_table_int(key, config_table_, "seq_block config");
}

double seq_block::retrieve_config_table_double(string key) const {
  return retrieve_table_double(key, config_table_, "seq_block config");
}

size_t seq_block::retrieve_config_table_size_t(string key) const {
  return retrieve_table_size_t(key, config_table_, "seq_block config");
}

string seq_block::name() const {
  return retrieve_config_table_str("name");
}

void seq_block::set_name(string new_name) {
  config_table_.set("name", new_name);
}

block_category seq_block::category() const {
  return retrieve_config_table("category").as<block_category>();
}
double seq_block::width_in_ms() const {
  return timeline2ms(timer_.width);
}

timeline seq_block::WIDTH() const {
  return timer_.width;
}

int seq_block::switch2loop(int index) {
  if (index == -1 || loop_ctrl_.loop_count == 1)
    return -1;

  int id = (index + 1) % loop_ctrl_.loop_count;
  if (id == 0)
    id = loop_ctrl_.loop_count;
  return id;
}

void seq_block::assign() {
}

void seq_block::local_loop_evolution() {
  if (loop_ctrl_.style == _local_array)
    for (int i = 0; i < loop_ctrl_.loop_count; i++)
      evolution(i);
}

vector<timeline> seq_block::absolute_keys() const {
  int n = timer_.keys.size();
  tlvec t = start_timeline() * tlvec::Ones(n);
  t += timer_.keys;
  return vector<timeline>(&t.data()[0], &t.data()[n]); // Note 0 ~ n-1.
}

bool seq_block::adjust_key_interval(timeline &t0, timeline &t1) const {
  t0 -= start_timeline();
  t1 -= start_timeline();
  if (t1 <= 0 || t0 >= WIDTH())
    return false; // outside of block width range.
  if (t0 < 0)
    t0 = 0;
  if (t1 > WIDTH())
    t1 = WIDTH();
  return true;
}

//bool seq_block::check_keys(timeline key0, timeline key1) const {
//    if (key0 >= 0 && key0 < key1 && key1 <= timer_.duration)
//        return true;
////  cout
////      << boost::format("%s %s %d %d.\n") % "S-S-L warning: "
////          % "invalid adjacent time keys (us): " % key0 % key1;
//    return false;
//}
void seq_block::write(ostream &ostr) const {
  ostr << get_header() << "\n";
}

string seq_block::get_header() const {
  string s;
  s += "# " + name() + "\n";
  s += "# " + to_string(loop_ctrl_.loop_prior) + "\n";
  s += "# width: " + to_string(width_in_ms()) + " ms.";
  return s;
}

void seq_block::plot() const {
}
} /* namespace seq */
} /* namespace ssl */
