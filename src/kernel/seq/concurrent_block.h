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

class concurrent_block :
    public glue {
 public:
  concurrent_block();
  ~concurrent_block();
  inline concurrent_block *Clone() const {
    return (new concurrent_block(*this));
  }
  virtual void evolution(int index = -1);
  virtual void write(std::ostream &ostr = std::cout) const;
  // should be used only prior to the seq{} formed.
  void sync();
  void allign();
 protected:
  virtual void assign();
 private:
  // used when one or more sub-blocks are concurrent block type.
  std::vector<seq_block *> simplify(const seq_block *sb);
};

}
}