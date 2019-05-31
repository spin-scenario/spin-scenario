Installation
=============
* **Binary installers**     

    Download the Spin-Scenario installers from the `releases`_. 

    .. note::
	    * To enable visualization, `Gnuplot <http://www.gnuplot.info/>`_ needs to be installed.
	    * To enable gpu acceleration, `ArrayFire <https://arrayfire.com/download/>`_ needs to be installed.
    
    * For Windows:     
    
        * run ``spin-scenario-1.0.0-win64.exe`` with the default installation option.  

        To run ``spin-scenario`` from any folder，you need to add ``C:\Program Files\spin-scenario-1.0.0\bin`` to system PATH.

    * For Ubuntu 16.04, 18.04:     
    
        .. code-block:: sh

         sudo dpkg -i spin-scenario-1.0.0-ubuntu.deb    
    
        Use ``sudo apt-get install -f`` to solve dependencies if the ``.deb`` package installation failed.
    
    * For CentOS 7:  
    
        .. code-block:: sh

         yum install boost-devel fftw-devel hdf5-devel NLopt-devel gnuplot		
         rpm -ivh spin-scenario-1.0.0-centos7.rpm --nodeps --force	
         ldconfig    	       

  
* **Build from source**

    * Clone source code:

        .. code-block:: sh

         git clone git@github.com:spin-scenario/spin-scenario.git


    * For Windows:

        .. note::

            Make sure you have ``visual studio`` and  ``cmake`` installed prior to the build. Download and install the precompiled 
            `Boost <https://sourceforge.net/projects/boost/files/boost-binaries/>`_ and `ArrayFire <https://arrayfire.com/download/>`_.  Then open ``CMakeLists.txt``:         
            
            * modify ``BOOST_ROOT`` into your Boost's full path.
            * modify ``ArrayFire_DIR`` into your ArrayFire's ``cmake`` directory (e.g. ``C:/Program Files/ArrayFire/v3/cmake``).

        * run ``sh ./dep_win10.sh`` to  generate visual studio projects for other 3rd-party libraries and ``spin-scenario``.   
        * complie 3rd-party libraries (``nlopt``, ``linenoise``, ``yacas``, ``lua``) via visual studio.   
        * complie ``spin-scenario`` via visual studio.   


    
    * For Linux:    

        Using the following commands to complie and install dependencies (Ubuntu).

        .. code-block:: sh

         sudo apt-get install libboost-all-dev libfftw3-dev libhdf5-dev libnlopt-dev libpython3-dev gnuplot
         sh ./dep.sh
        
        Using the following commands to complie and install spin-scenario.   

        .. code-block:: sh

         mkdir build
         cd build
         cmake ..
         make
         make install



.. _releases: https://github.com/spin-scenario/spin-scenario/releases

