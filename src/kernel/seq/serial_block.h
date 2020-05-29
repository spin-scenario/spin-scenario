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
#include "glue.h"

namespace ssl {
namespace seq {

class serial_block :
    public glue {
 public:
  serial_block();
  ~serial_block();
  inline serial_block *Clone() const {
    return (new serial_block(*this));
  }
  virtual void evolution(int index = -1);
  virtual void write(std::ostream &ostr = std::cout) const;
 protected:
  virtual void assign();
 private:
  std::vector<const seq_block *> retrieve_cycle_seq_blocks(int prior) const;
  int max_cycle_count(int prior) const;
  imat descending_cycle_matrix() const;
  std::vector<std::vector<int>> cycle_index_list(const std::vector<std::vector<int>> &old, int count) const;
};

}
}