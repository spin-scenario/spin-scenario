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

-- example: spectrum of two spin system (strong coupling). 

-- scenario A: spin system generation.
B0{"500 MHz"}
local sys =spin_system{
  spin = '1H 1H',
    zeeman = '1 scalar 0.5 ppm 2 scalar 1 ppm',
    jcoupling = '1 2 scalar 7 Hz'
}

-- scenario B: pulse sequence assembly.
local rf45 =hardRF{beta =45}
local adc =acq{np =1024, sw =10000}

local fid =seq{rf45, adc}
-- scenario C: experimental study.
result =run{exp =fid, spinsys =sys}


local t = linspace(0,1024/10000*1e3,1024) -- acquisition time axis, unit ms.
local x = linspace(-5000,5000,1024)    -- spectrum axis, unit Hz.

plot("title<fid:re> xlabel<time/ms>",lines(t, result["fid:re"]))  
plot("title<spec:abs> xlabel<freq/Hz> xrange<-1000:1000>",lines(x, result["spec:abs"]))
