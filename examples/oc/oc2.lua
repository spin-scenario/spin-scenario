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

-- example: broadband pulse optimization for proton.
B0{"500 MHz"}
-- set up spin system.
local config = {
    spin = '1H',
    zeeman = '1 scalar -5:5:41 kHz' -- this will create 41 chemical shift values within [-5,5] kHz.
}

local sys = spin_system(config)
local oc = rf_optimizer(sys)
local rf90 = oc:optimize{ width = 5.12, step = 512, init_pattern = 'rand', init_state = 'I1z', targ_state = 'I1x' }

plot(rf90)
--rf90:switch('ux/uy')
--plot(rf90)

oc:projection{init_state = 'I1z', rf = rf90, observ_states = {'I1z', 'I1x', 'I1y'}, option ='broadband'}

--specgram(rf90, { window = 'hamming', wlen = 128, overlap = 0.98, nfft = 1024, style = '' })
