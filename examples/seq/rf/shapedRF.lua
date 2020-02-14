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

-- example: shaped pulses. 

local rf1 = shapedRF{width =5.12, step =256, max_amp =100, pattern ="sinc", lobe =7}
write("sinc.RF",rf1)

-- by default, the RF data file contains two columns in amp(Hz)/phase(deg).
-- Note if the raw shape is in ux(Hz)/uy(Hz), please use additional option mode = 'ux/uy'.
local rf2 = shapedRF{width = 5.12, pattern = "shape.RF"} 

local rf3 = shapedRF{width =10, step =100, pattern ="2*(10-t) + 5*Cos(10*t)^2"} 
-- turn to YACAS for more usage to write your expression.

plot(rf1, rf2, rf3) 