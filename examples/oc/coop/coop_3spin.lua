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

-- example: cooperative pulses optimization (three-qubit pseudo-pure state preparation).
--  Appl. Phys. Lett. 104, 242409 (2014); https://doi.org/10.1063/1.4884295 

B0{"500 MHz"}
-- set up spin system.
local config = {
    spin = '13C 13C 13C',
    zeeman = '1 scalar 15.74 kHz 3 scalar -4.302 kHz',
    jcoupling = '1 2 scalar 54.2 Hz 2 3 scalar 35.4 Hz'
}
local sys = spin_system(config)

local oc = multi_rf_optimizer(sys)
local opt_pulse = oc:optimize {
    ncoop = 7,
    width = 18.5,
    step = 1850,
    init_pattern = 'rand',
    init_state = 'I1z+I2z+I3z',
    targ_state = 'I1z+I2z+I3z+2*I1zI2z+2*I1zI3z+2*I2zI3z+4*I1zI2zI3z',
    tol_f = 1e-9,
    max_eval = 200
}

oc:projection{init_state = 'I1z+I2z+I3z', coop_rf = opt_pulse, observ_states = { 'I1z', 'I2z', 'I3z', '2*I1zI2z', '2*I1zI3z', '2*I2zI3z', '4*I1zI2zI3z' }}

for i = 1, #opt_pulse do
    opt_pulse[i]:switch('ux/uy')
    plot(opt_pulse[i])
    --specgram(opt_pulse[i], { wlen = 128, overlap = 0.96, nfft = 1024, style = '' })
end
