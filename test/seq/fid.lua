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

-- example: FID sequence. 

-- scenario A: spin system generation.
B0{"500 MHz"}
local proton =spin_system{
    spin ="1H",
    zeeman ="1 scalar 100 Hz",
    relaxation ="T1 1000 ms T2 100 ms"
}
-- scenario B: pulse sequence assembly.
local rf45 =hardRF{beta =45}
local adc =acq{np =4096, sw =10000, phase= "-y"}

local fid =seq{rf45, adc}
-- scenario C: experimental study.
result =run{exp =fid, spinsys =proton} -- the var result is global so that we can still deal with it on the terminal when seq simulation done. 


local t = linspace(0,4096/10000*1e3,4096) -- acquisition time axis, unit ms.
local x = linspace(-5000,5000,4096)    -- spectrum axis, unit Hz.

plot("title<fid:re> xlabel<time/ms> gnuplot<unset key>",lines(t, result["fid:re"]))  
plot("title<spec:abs> xlabel<freq/Hz> xrange<-500:500> gnuplot<unset key>",lines(x, result["spec:abs"]))

