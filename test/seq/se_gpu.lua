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

-- example: spin echo imaging. 

-- basic parameters for imaging.
B0{"3 T"}
peak_grad{40} -- T/m
slew_rate{200} -- T/m/s
pw90{5} -- us
seq_parm{fov = '240*240', matrix = '128*128'}
reduce_phantom{z0=180, z1=185, dx=2, dy=2} -- SEE NOTE for detail. 
-- The above line will reduce the spin ensemble number to those slices z0~z1, and for each slice spins every dx and every dy were used for simulation.
-- Valid options include x0, x1, y0, y1, dx, dy, dz. In this way the simulations will be faster, but may lead to poor image quality. 
-- Commit this if you wish to get better image results, which may lead to longer computing time.
-- This API will be removed in future.  

-- pulse sequence assembly.
local gy =trapGrad{axis ="Y", func ="phase_encode", width =2}
local adc =acq{np =128, sw =32000}
local gx =trapGrad{axis ="X", func ="read_out"}
local gxPre =trapGrad{axis ="X", area =0.5*area(gx), width =2}

local gyspoil =trapGrad{axis ="Y", area = math.pi*2*1e3/42.57/(240/128), width =3}

local rf90 =hardRF{beta=90, width=0.002}
local rf180 =hardRF{beta=180, width=0.002}
local TR =500
local TE = 15

d1 = delay{width = TE/2-rf90.tau/2-gy.tau-rf180.tau/2}
d2 = delay{width = TE/2-gx.tau/2-rf180.tau/2}
d3 =delay{width= TR-TE-gx.tau/2-rf90.tau/2-gyspoil.tau}

local se =seq{rf90, d1, gxPre + gy#, rf180, d2, gx + adc, d3, gyspoil}

-- experimental study.
result =run{exp =se, phantom ="mida.h5", acc="gpu"} 

-- the 128*128 2D FFT data was saved into matrix such as result["IMG:re"] or result["IMG:im"] or result["IMG:abs"].
plot('color<Greys>',map(result["IMG:abs"]))
write('img_se.txt',result["IMG:abs"])


-- NOTE. 
-- This is a temporary API used to reduce the number of total spins involved in the simulation in order to perform a quick test. Commit this line will use all spins inside the phantom. ONLY spins in silces [z0-z1] are involved in the imaging.
-- MIDA phantom: X/Y/Z 480*480*350, spatial resolution: 0.5 mm 
--  MNI phantom: X/Y/Z 90*108*90, spatial resolution: 2 mm 
-- The three-dimensional origin locates at the center of phantoms by default.

