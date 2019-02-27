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

-- example: pulse optimization for two homonuclear spin system.
B0{"500 MHz"}
-- set up spin system.
local config = {
    spin = '1H 1H',
    zeeman = '1 scalar 0.898 ppm',
    jcoupling = '1 2 scalar 10.2 Hz'
}
local sys = spin_system(config)

local oc = rf_optimizer(sys)
local opt_pulse = oc:optimize{width = 24.5, step = 512, init_state = 'I1xI2x', targ_state = 'I1yI2y', max_eval = 50}

plot(opt_pulse)
