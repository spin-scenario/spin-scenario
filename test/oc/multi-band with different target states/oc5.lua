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

-- example: pulse optimization for multi-band with different target states.  ---````--- (e.g. side band, center band, side band)
B0{"500 MHz"}
local config = {
    spin = '1H',
 zeeman = '1 scalar -1:1:288 kHz' -- this will create 288 chemical shift values within [-1,1] kHz.
}
local sys = spin_system(config)

-- define targ state profile.
local Ix, Iy, Iz = sys:expr_state('I1x'), sys:expr_state('I1y'), sys:expr_state('I1z')

-- [-1,1]kHz [central band (1/3 width) suppressed and side bands (1/3 width) excited.]
local np = 100 -- discrete points for each band.
local targ_list = {}

for i = 1, np - 1 do targ_list[i] = Ix end 
for i = 1, np do targ_list[i + np - 1] = Iz end  -- central band
for i = 1, np - 1 do targ_list[i + 2 * np - 1] = Ix end
-- total 288 points match with spin sys define. 3np-2=288


local oc = rf_optimizer(sys)
local opt = oc:optimize{ width = 10, step = 500, init_pattern = 'rand', init_state = 'I1z', targ_state = targ_list }
--print('average power: ' .. opt.average_power)
--os.execute("mkdir cacl")
--ssl.write('cacl/rf_shape.txt', opt)
--opt:switch('amp/phase')
plot(opt)
oc:projection { init_state = 'I1z', rf = rf90, observ_states = {}, option = 'broadband'}
