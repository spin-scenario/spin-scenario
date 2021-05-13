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

-- example: slice view of phantom model.
--local mida = phantom("mida.h5", "mida_config.txt")
--mida:view{x=240, y=240, z=175} 
--mida:view{prop="T1", z=175} 
-- MIDA phantom: X/Y/Z 480*480*350, spatial resolution: 0.5 mm 

--local p1 = phantom("circles.h5", "config.txt") --a 2d phantom created by h5_phantom_2d.m
--p1:view{prop="T1", z=1} 

--local p2 = phantom("spheres.h5", "config.txt") --a 3d phantom created by h5_phantom_3d.m
--p2:view{x=100, y=30, z=50} 