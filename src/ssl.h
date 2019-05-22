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
#include <kernel/spinsys/spin_system.h>
using namespace ssl::spinsys;
#include <kernel/sample/phantom.h>
using namespace ssl::sample;
#include <kernel/seq/pulse_seq.h>
#include <kernel/seq/seq_block_factory.h>
using namespace ssl::seq;
#include <kernel/oc/grape.h>
#include <kernel/oc/coop_grape.h>
#include <kernel/oc/tf_opt.h>
using namespace ssl::oc;
#include <kernel/utilities/ssl_plot.h>
using namespace ssl::utility;

namespace ssl {

void bindings(sol::state &lua);

}
