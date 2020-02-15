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

-- example: vlidation of the result from "coop_2spin.lua" (two-qubit pseudo-pure state preparation).
output_terminal{type = "eps", font = "Arial,16"} -- available font includes png, eps and tex.
B0{"500 MHz"}
-- set up spin system.
local config = {
    spin = '13C 1H',
    jcoupling = '1 2 scalar 140 Hz'
}
local sys = spin_system(config)



local opt_pulse = multiRF{channel = '13C 1H', width = 3.57, step = 714, file = 'coop_20190221_113357.h5', mode = 'amp/phase' }

local oc = multi_rf_optimizer(sys)
oc:projection{init_state = 'I1z+4*I2z', coop_rf = opt_pulse, observ_states = { 'I1z', 'I2z', '2*I1zI2z' } }

for i = 1, #opt_pulse do
    opt_pulse[i]:switch('ux/uy')
    plot(opt_pulse[i])
    --specgram(opt_pulse[i], { wlen = 64, overlap = 0.95, nfft = 2048, style = '' })
end
