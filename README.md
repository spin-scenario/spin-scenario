<a><img src="docs/source/media/logo.png" width="400"></a>

Features
-------------------------------
The Spin-Scenario is an intuitive, flexible and unique scripting environment for realistic MR simulations.     

* **General, fast MR simulations**   
    It uses accelerated Liouville space computing model being compatible with MR imaging and MR spectroscopy.

* **Flexible scenario scripting**    
    The powerful Lua bindings offer users a flexible, intuitive and unique scripting environment for the creation of MR scenarios.

* **Elegant pulse sequence programming**   
    The specially designed programming syntax enables a clear, concise, and accurate description of pulse sequences.  

* **Efficient pulse optimization**   
    The build-in optimal control module provides an easy way for optimization of shaped pulses as well as cooperative pulses.

Quick Start
--------------------------------------

  * Open the terminal and start the environment with `spin-scenario`.
  * To run a scenario script, simply use the command like `load('fid.lua')`. 
  * To switch to old script, use the up and down key.  
  * To quit the program, use `q`.
  
<p align="center">
<img src="docs/source/media/cmd_screen.png"  width="480">
</p>

Alternatively, users may also execute the scenario script in terminal like ``spin-scenario fid.lua``, just make sure that the script file locates at the current path of the terminal.

As a start, we show a simple [FID](examples/seq/acrylic.lua) sequence as follow:  
```lua
        -- scenario A: spin system generation.
        B0{"500 MHz"}
       local acrylic =spin_system{
                spin = "1H 1H 1H",
              zeeman = "2 scalar 88.42 Hz 3 scalar 214.9 Hz",
           jcoupling = "1 2 scalar 10.4 Hz 1 3 scalar 1.2 Hz 2 3 scalar 17.4 Hz"}
        -- scenario B: pulse sequence assembly.
        local rf45 =hardRF{beta =45}
        local adc =acq{np =1024, sw =500}

        local fid =seq{rf45, adc}
        -- scenario C: experimental study.
        result =run{exp =fid, spinsys =acrylic}
```
<a><img src="docs/source/media/seq_fid_acrylic_signal.png" width="360"></a> <a><img src="docs/source/media/seq_fid_acrylic_spec.png" width="360"></a> 

More scenario examples can be found in [examples](examples). 

Documentation
--------------------------------------
Spin-Scenario Documentation [http://spin-scenario.rtfd.io/](http://spin-scenario.rtfd.io/).

Contributing to Spin-Scenario
---------------------------------------
Spin-Scenario is provided as an open-source project under the `Apache License 2.0`, and will be continuously developed for highly usability as well as high-performance.
    
Spin-Scenario is mainly written and maintained by Yan CHANG (`changy@sibet.ac.cn`). To make it better, you are welcome to contribute in many ways to help: helping with bugs and feature requests, discussing the design and the API, uploading scenario scripts and writing low-level code.
    
   To cite Spin-Scenario using a BibTeX entry:
                    
        @article{CHANG2019,
        title = "Spin-Scenario: A flexible scripting environment for realistic MR simulations",
        journal = "Journal of Magnetic Resonance",
        volume = "301",
        pages = "1 - 9",
        year = "2019",
        issn = "1090-7807",
        doi = "https://doi.org/10.1016/j.jmr.2019.01.016",
        url = "http://www.sciencedirect.com/science/article/pii/S1090780719300229",
        author = "Yan Chang and Daxiu Wei and Huihui Jia and Cecilia Curreli and Zhenzhou Wu and Mao Sheng and Steffen J. Glaser and Xiaodong Yang"
        }

 <p align="center">
<img src="docs/source/media/cover.png"  width="480">
</p>

 Acknowledgments
---------------------------------------
This work was supported by National Natural Science Foundation of China 11505281, 11675254.
