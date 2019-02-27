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

-- example: a modified INEPT arrayed experiments where the S-spin signal is observed immediately after it has been transferred from I for a set of varied time in delay block tau. 

-- scenario A: spin system generation.
B0{"500 MHz"}

local sys = spin_system {
    spin = '1H 13C',
    zeeman = '1 scalar -100 Hz 2 scalar 500 Hz',
    jcoupling = '1 2 scalar 140 Hz'
}

-- scenario B: pulse sequence assembly.
local rf90I = hardRF{ channel = '1H', beta = 90, phase = 'x' }
local rf180 = hardRF{ channel = '1H|13C', beta = 180, phase = 'x|x' }
local  rf90 = hardRF{ channel = '1H|13C', beta = 90, phase = 'y|y' }

local tau = delay { width = '1.7857:3:8' } -- note here the 'delay' block was varied 8 times.
local adc = acq{np = 1024, sw = 10000, channel ='13C'}

local inept = seq{rf90I, tau#, rf180, tau# , adc}

-- scenario C: experimental study.
result =run{exp =inept, spinsys =sys}


local t = linspace(0,1024/10000*1e3,1024) -- acquisition time axis, unit ms.
local x = linspace(-5000,5000,1024)    -- spectrum axis, unit Hz.

plot("title<fid:re> xlabel<time/ms>",lines(t, result["fid:re"]))  
plot("title<spec:abs> xlabel<freq/Hz> xrange<0:1000>",lines(x, result["spec:abs"]))

