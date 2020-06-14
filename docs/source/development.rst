Development
=============
  
* **Build from source**

    * Clone source code:

        .. code-block:: sh

         git clone git@github.com:spin-scenario/spin-scenario.git


    * For Windows:

        .. note::

            Make sure you have `Visual Studio (2019) <https://visualstudio.microsoft.com>`_ and `Cmake <https://cmake.org/download/>`_ installed prior to the build. 
            
            Download and install the precompiled 
            `Boost (1.7.2) <https://sourceforge.net/projects/boost/files/boost-binaries/>`_, `CUDA <https://developer.nvidia.com/cuda-downloads>`_, `ArrayFire <https://arrayfire.com/download/>`_
            and `HDF5 <https://www.hdfgroup.org/downloads/hdf5>`_.  
            
            Then open ``CMakeLists.txt``:
                       
            * modify ``BOOST_ROOT`` into your Boost's full path.
            * modify ``ArrayFire_DIR`` into your ArrayFire's ``cmake`` directory (e.g. ``C:/Program Files/ArrayFire/v3/cmake``).

        * run ``sh ./dep_win.sh`` to  generate visual studio projects for other 3rd-party libraries (``nlopt``, ``linenoise``, ``yacas``, ``lua``) and ``spin-scenario``.   
        * complie the above 3rd-party libraries via visual studio.   
        * complie ``spin-scenario`` via visual studio.   


    
    * For Linux:    

        Using the following commands to complie and install dependencies (Ubuntu).

        .. code-block:: sh

         sudo apt-get install cmake libboost-all-dev libfftw3-dev libhdf5-dev libnlopt-dev libnlopt-cxx-dev libpython3-dev gnuplot
         sh ./dep.sh
        
        Using the following commands to complie and install spin-scenario.   

        .. code-block:: sh

         mkdir build
         cd build
         cmake ..
         make
         make package
         make install



.. _releases: https://github.com/spin-scenario/spin-scenario/releases


* **Documentation**     

    Spin-Scenario use `Sphinx <http://www.sphinx-doc.org/en/master/index.html>`_ to create all the documentation.
    Go into the `docs`_ directory and use ``sphinx-build -b html source build`` to build HTML files.

.. _docs: https://github.com/spin-scenario/spin-scenario/tree/master/docs
