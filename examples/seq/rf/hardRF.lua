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

-- example: Hard RF pulse. 

-- 60 degree -y-pulse for proton 1H excitation.
local rf1 = hardRF{beta =60, phase="-y"}
write("rf1.RF",rf1)

-- a refocusing pulse applied on a heteronuclear 1H-13C spin system
local rf2 = hardRF{beta =180, channel ="1H|13C", phase="x|x"}
write("rf2.RF",rf2)
