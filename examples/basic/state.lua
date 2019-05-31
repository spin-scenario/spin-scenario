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

-- example: show how to create the density state (similar for operator). 
B0{"500 MHz"}

-- set up spin system.
local config = {
    spin = '1H 1H 1H',
    zeeman = '2 scalar -38.8 Hz 3 scalar 412 Hz',
    jcoupling = '1 3 scalar 17 Hz 2 3 scalar 11 Hz'
}

local sys = spin_system(config)

local I1x = 0.5*(sys:state{ 1, 'I+' } + sys:state{ 1, 'I-' })
local I1x_ = sys:expr_state('I1x')
print(I1x_-I1x)  -- should result zero vector.


local I1xI3x = 0.25*(sys:state{ 1, 'I+', 3, 'I+' } + sys:state{ 1, 'I+', 3, 'I-' } + sys:state{ 1, 'I-', 3, 'I+' } + sys:state{ 1, 'I-', 3, 'I-'})
local I1xI3x_ = sys:expr_state('I1xI3x')

local I1yI3y = -0.25*(sys:state{ 1, 'I+', 3, 'I+'} - sys:state{ 1, 'I+', 3, 'I-'} - sys:state{ 1, 'I-', 3, 'I+'} + sys:state{ 1, 'I-', 3, 'I-'})
local I1yI3y_ = sys:expr_state('I1yI3y')

--print(I1xI3x-I1xI3x_)  -- should result zero vector.
--print(I1yI3y-I1yI3y_)  -- should result zero vector.

local rho = sys:expr_state('2*I1xI3x-0.5*I2xI3y')
print(rho)
write('rho.txt', rho)


