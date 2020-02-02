.. sol documentation master file, created by
   sphinx-quickstart on Mon Feb 29 21:49:51 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: media/spin-scenario.png
	:target: https://github.com/spin-scenario/spin-scenario
	:alt: spin-scenario repository
	:align: center

Spin-Scenario |version|
=======================
*a flexible, unified scripting environment for realistic MR simulations*

Get going:
----------

.. toctree::
	:maxdepth: 1
	:name: mastertoc
	
	about
	installation
	tutorial/tutorial-top
	module/spinsys
	module/phantom
	module/seq
	module/oc
	module/physx
	module/utility
	api/api-top
	development
	roadmap
	publication
	license
	contact

the basics:
-----------
.. code-block:: lua
	:linenos:

	 -- scenario A: spin system generation.
        B0{"500 MHz"}
        local proton =spin_system{
            spin ="1H",
            zeeman ="1 scalar 100 Hz",
            relaxation ="T1 1000 ms T2 100 ms"
        }
        -- scenario B: pulse sequence assembly.
        local rf45 =hardRF{beta =45}
        local adc =acq{np =4096, sw =10000}

        local fid =seq{rf45, adc}
        -- scenario C: experimental study.
        result =run{exp =fid, spinsys =proton}

|fid_raw| |fid_spec|

.. note::
	The code above and more examples can be found in the `examples directory`_.

Acknowledgments
---------------
Spin-Scenario was supported by National Natural Science Foundation of China 11505281, 11675254.

.. _spin-scenario: https://github.com/spin-scenario/spin-scenario
.. _issues: https://github.com/spin-scenario/spin-scenario/issues
.. _examples directory: https://github.com/spin-scenario/spin-scenario/tree/master/examples

.. |fid_raw| image:: media/seq_fid_raw.png
	:height: 320
	:align: middle

.. |fid_spec| image:: media/seq_fid_spec.png
	:height: 320
	:align: middle

.. |nsfc| image:: media/logo-NSFC.jpg
	:height: 45
	:align: middle	

.. |sibet| image:: media/sibet.png
	:height: 45
	:align: middle

.. |ecnu| image:: media/ecnu.gif
	:height: 40
	:align: middle	

.. |children_hosptial| image:: media/children_hosptial.png
	:height: 45
	:align: middle	

.. |tum| image:: media/tum.gif
	:height: 45
	:align: middle	
