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

-- example: gradient pulses. 


peak_grad{40} -- T/m
slew_rate{200} -- T/m/s
local grad = exprGrad{axis ='Y', expr ='Sin(t*7)', width =2}
write("grad.txt", grad)
plot(grad)

os.exit()
local g1 = trapGrad{axis ='x', width = 2, area = 70} 
local g2 = trapGrad{axis ='z', width = 2, amp = 35}
local g3 = trapGrad{axis ='x', flat_time = 5, amp = 40}
local g4 = trapGrad{axis ='y', flat_time = 5, flat_area = 150}

plot(g1, g2, g3, g4)

seq_parm{fov = '100*100', matrix = '32*32'}
local gy =trapGrad{axis ="Z", func ="phase_encode", width = 0.5}
plot(gy)

--local grad = exprGrad{axis ='Y', expr ='Sin(t*7)', width =2}



