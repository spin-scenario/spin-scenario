# Spin-Scenario

![Image text](https://github.com/spin-scenario/spin-scenario-doc/blob/master/source/graphical_abstract.png)
  
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

Installation on Linux
--------------------------------------
* **Binary installers**     

    Download the Spin-Scenario installers from the [releases](https://github.com/spin-scenario/spin-scenario/releases).  
    * For Ubuntu 16.04, 18.04:     
    
        * `sudo dpkg -i spin-scenario-1.0.0-ubuntu16.04.deb`     
        * `sudo dpkg -i spin-scenario-1.0.0-ubuntu18.04.deb`                         
    
    Use `sudo apt-get install -f` to solve dependencies if the `.deb` package installation failed.
    
    * For CentOS 7:  
    
         * `yum install boost-devel fftw-devel hdf5-devel NLopt-devel gnuplot`		
         * `rpm -ivh spin-scenario-1.0.0-centos7.rpm --nodeps --force`	
         * `ldconfig`       	       

  
  
* **Build from source**

    The Spin-Scenario was developed based on a list of opensource libraries, which need to be complied in advanced. 
    To do this, first install some packages by `sudo apt-get install libboost-all-dev libfftw3-dev libhdf5-dev libnlopt-dev libpython3-dev` (Ubuntu). 
    Then run `sh ./dep.sh` in the root of the repository to complie other packages, and this should work automatically.
    
    Using the following commands to complie and install spin-scenario.
      
    `mkdir build`   
    `cd build`  
    `cmake ..`   
    `make`  
    `make install`
  
Basic usage
--------------------------------------

  * Open the terminal and start the environment with `spin-scenario`.
  * To run a scenario script, simply use the command `load('demo.lua')`. 
  * You may also run script snippet on the terminal, e.g. to process the results manually after the scenario playing. 
  * To switch to old script, just the up and down key.  
  * To quit the program, use `q`.	    
  * Scenario examples can be found in [test](https://github.com/spin-scenario/spin-scenario/tree/master/test).
  
       
        -- A FID sequence demo for acrylic spin system.
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
        result =run{exp =fid, spinsys =acrylic}
   
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
        
![Image text](https://github.com/spin-scenario/spin-scenario-doc/blob/master/source/cover.png)

 Acknowledgments
---------------------------------------
This work was supported by National Natural Science Foundation of China 11505281, 11675254.
