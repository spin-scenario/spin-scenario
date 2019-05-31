--[[ Copyright 2019 The Spin-Scenario Authors. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================-]]

-- preliminary test for tensorflow oc of one spin system. 

-- NOTE: to test this script, you need to install tensorflow into the system. 
-- sudo apt install python3-dev python3-pip
-- sudo -H pip3 install tensorflow==1.10
-- see https://tensorflow.google.cn/install/pip

B0{"500 MHz"}
-- set up spin system.
local config = {
    spin = '1H'
    --zeeman = '1 scalar -2.898 ppm'
}
local sys = spin_system(config)
local oc = rf_optimizer_tf(sys)
local opt_pulse = oc:optimize{ width = 1, step = 100, init_pattern = 'rand', init_state = 'I1z', targ_state = '-I1x'}

plot(opt_pulse)

oc:projection{init_state = 'I1z', rf = opt_pulse, observ_states = {'I1z', 'I1x', 'I1y'} }
