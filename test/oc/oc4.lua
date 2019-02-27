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

-- example: pulse optimization for two heteronuclear spin system.
B0{"500 MHz"}
-- set up spin system.
local config = {
    spin = '1H 13C',
    jcoupling = '1 2 scalar 140 Hz',
}
local sys = spin_system(config)

local oc = rf_optimizer(sys)
local opt_pulse = oc:optimize {
    width = 7.14,
    step = 357,
    init_pattern = 'rand',
    init_state = 'I1x',
    targ_state = 'I2x',
    --tol_f = 1e-8,
    max_eval = 100
}

opt_pulse:set_name('I_x to S_x')
opt_pulse:switch('ux/uy')
plot(opt_pulse)

oc:projection{ init_state = 'I1x', rf = opt_pulse, observ_states = {} }
