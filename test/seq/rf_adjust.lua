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

-- example: shaped RF scaled to specific flip angle (180 degree). 

-- scenario A: spin system generation.
B0{"500 MHz"}
local proton =spin_system{
spin ="1H",
relaxation="T1 1000 ms T2 20 ms"
}
-- scenario B: pulse sequence assembly.
local rf180 =shapedRF{pattern="sinc", width= 2, step=100, gain="10:2500:100"}
plot(rf180) -- show the standard sinc pulse above.

local adc =acq{np =1024, sw =10000}
local d1 = delay{width = 1000}

local fid =seq{rf180#, adc, d1}

-- scenario C: experimental study.
result =run{exp =fid, spinsys =proton}

-- each row of the fid data was stored into tables such as result["fid:re"] or result["fid:im"] or result["fid:abs"], for this demo, 50 acquisitions were taken, thus each table contains 50 raw vectors.

-- the corresponding 1d spec of the above raw data was stored into tables such as result["spec:re"] or result["spec:im"] or result["spec:abs"].

-- all the above 1d specs were stored into a matrix of which each row corresponds to a specific rf180 amplitude gain. The matrix can be obtained by
-- result["SPEC:re"] or result["SPEC:im"] or result["SPEC:abs"]
-- using the result["SPEC:abs"], it is easier to find the first trough in the SPEC, which is the gain for 180 degree sinc pulse.

plot("title<1d spec row-by-row> xlabel<time/ms> ylabel<gain value>",map(result["SPEC:abs"], "style<3d> xrange<0:102.4> yrange<10:2500>"))
-- use style<image> to show the map plot of the spec. 
