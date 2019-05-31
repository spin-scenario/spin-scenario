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

local mida = phantom('mida.h5')
mida:view{z = 150, x= 120, y = 180} -- this will show the T1 map for slice No. 150 in Z axis, No. 120 in X axis and No. 180 in Y axis respectively.


--local mni = phantom('mni.h5')
--mni:view{z = 45, x= 50, y = 60}


-- MIDA phantom: X/Y/Z 480*480*350, spatial resolution: 0.5 mm 
--  MNI phantom: X/Y/Z 90*108*90, spatial resolution: 2 mm 
