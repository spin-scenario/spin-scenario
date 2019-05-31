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

-- example: spin system operations.
B0{"500 MHz"}
-- set up spin system.
local sys = spin_system{
    spin = '1H 13C',
    jcoupling = '1 2 scalar 140 Hz',
    relaxation = 'T1 1350 1000 ms T2 200 50 ms'
}

local Ix = sys:expr_state("I1x")
local Sx = sys:expr_state("I2x")
local IxpSx = sys:expr_state("I1x+I2x")
print(Ix,Sx)
write('s1.txt', Ix, Sx, IxpSx)


local R = sys.R
local L0 = sys.L0
local L = sys.L
write('s2.txt', R, L0, L)





