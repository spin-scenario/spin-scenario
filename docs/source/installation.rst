Installation
=============
   
    .. note::
	    To install from source, refer to :doc:`Development<development>` for details.

* **Binary installers**     

    Download the Spin-Scenario installers from the `releases`_. 

    .. note::
	    * To enable visualization, `Gnuplot <http://www.gnuplot.info/>`_ needs to be installed.
	    * To enable gpu acceleration, `ArrayFire <https://arrayfire.com/download/>`_  and `CUDA <https://developer.nvidia.com/cuda-downloads>`_ need to be installed.
    
    * For Windows:     
    
        * run ``Spin-Scenario-1.1.0-win64.exe`` to start installation.  

        To run ``spin-scenario`` from any folder，you need to add your ``spin-scenario`` binary path e.g. ``C:\Program Files\Spin-Scenario\bin`` to system PATH environment.

        .. note::
            You may need to install `Visual C++ Redistributable Package <https://aka.ms/vs/16/release/vc_redist.x64.exe>`_ in case of encountering missing dll file such as **vcruntime140_1.dll** when start the Spin-Scenario.

    * For Ubuntu 16.04, 18.04, 20.04:     
    
        .. code-block:: sh

         sudo dpkg -i Spin-Scenario-1.1.0-Ubuntu.deb    
    
        Use ``sudo apt-get install -f`` to solve dependencies if the ``.deb`` package installation failed.
    
    * For CentOS 7:  
    
        .. code-block:: sh

         yum install boost-devel fftw-devel hdf5-devel NLopt-devel gnuplot		
         rpm -ivh spin-scenario-1.0.0-centos7.rpm --nodeps --force	
         ldconfig    	       
  

.. _releases: https://github.com/spin-scenario/spin-scenario/releases

