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

-- example: spectrum of acrylic acid. 

-- scenario A: spin system generation.
B0{"500 MHz"}
local acrylic =spin_system{
    spin ="1H 1H 1H",
    zeeman ="2 scalar 88.42 Hz 3 scalar 214.9 Hz",
    jcoupling ="1 2 scalar 10.4 Hz 1 3 scalar 1.2 Hz 2 3 scalar 17.4 Hz"
}
-- scenario B: pulse sequence assembly.
local rf45 =hardRF{beta =45}
local adc =acq{np =1024, sw =500}

local fid =seq{rf45, adc}
-- scenario C: experimental study.
result =run{exp =fid, spinsys =acrylic} -- the var result is global so that we can still deal with it on the terminal when seq simulation done. 


local t = linspace(0,1024/500*1e3,1024) -- acquisition time axis, unit ms.
local x = linspace(-250,250,1024)    -- spectrum axis, unit Hz.

plot("title<fid:re> xlabel<time/ms>",lines(t, result["fid:re"]))  
plot("title<spec:abs> xlabel<freq/Hz> xrange<-250:250>",lines(x, result["spec:abs"]))

-- each row of the fid data was stored into tables such as result["fid:re"] or result["fid:im"] or result["fid:abs"], for this demo, only one acquisition was taken.

-- the corresponding 1d spec of the above raw data was stored into tables such as result["spec:re"] or result["spec:im"] or result["spec:abs"].

