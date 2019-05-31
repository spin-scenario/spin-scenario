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

-- example: COSY test. 

-- scenario A: spin system generation.
B0{"500 MHz"}

local sys =spin_system{
    spin ="1H 1H 1H",
    zeeman ="2 scalar 1000 Hz 3 scalar -1000 Hz",
    jcoupling ="1 2 scalar 200 Hz 1 3 scalar 100 Hz 2 3 scalar 150 Hz",
    relaxation ="T1 1000 1200 900 ms T2 20 15 35 ms"
}

-- scenario B: pulse sequence assembly.
local rf90 =hardRF{beta=90, phase="x"}
local rf180 =hardRF{beta=180, phase="y"}
local d1=delay{width="0:31.75:128"}  -- note the increment delay equals 1/sw
local d2=delay{width=1000}
local adc =acq{np =128, sw =4000}

local cosy =seq{rf90, d1#, rf180, adc, d2}

-- scenario C: experimental study.
result =run{exp =cosy, spinsys =sys}  -- the var result is global so that we can still deal with it on the terminal when seq simulation done. 


-- each row of the fid data was stored into tables such as result["fid:re"] or result["fid:im"] or result["fid:abs"], for this demo, 128 acquisitions were taken, thus each table contains 128 raw vectors.

-- the corresponding 1d spec of the above raw data was stored into tables such as result["spec:re"] or result["spec:im"] or result["spec:abs"].

-- the 2D FFT data was saved into matrix such as result["IMG:re"] or result["IMG:im"] or result["IMG:abs"].

plot("title<2D spec> xlabel<F1/Hz> ylabel<F2/Hz>",map(result["IMG:abs"], "style<image> xrange<-2000:2000> yrange<-2000:2000>"))
plot("title<2D spec> xlabel<F1/Hz> ylabel<F2/Hz>",map(result["IMG:abs"], "style<3d> xrange<-2000:2000> yrange<-2000:2000>"))
-- use style<3d> to show the surface plot of the spec. 


