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

-- example: cooperative pulses optimization (two-qubit pseudo-pure state preparation).

B0{"500 MHz"}
-- set up spin system.
local config = {
    spin = '13C 1H',
    jcoupling = '1 2 scalar 140 Hz'
}
local sys = spin_system(config)

local oc = multi_rf_optimizer(sys)
local opt_pulse = oc:optimize{
    ncoop = 3,
    width = 3.57,
    step = 714,
    init_pattern = 'rand',
    init_state = 'I1z+4*I2z',
    targ_state = 'I1z+I2z+2*I1zI2z',
    tol_f = 1e-9,
    max_eval = 100
}

oc:projection{init_state = 'I1z+4*I2z', coop_rf = opt_pulse, observ_states = { 'I1z', 'I2z', '2*I1zI2z'}}

for i = 1, #opt_pulse do
    opt_pulse[i]:switch('ux/uy')
    plot(opt_pulse[i])
    --specgram(opt_pulse[i], { wlen = 64, overlap = 0.95, nfft = 2048, style = '' })
end
