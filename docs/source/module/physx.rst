
************
Physx Engine
************

When both :doc:`pulse sequence<seq>` as well as the :doc:`spin system<spinsys>` (or :doc:`phantom<phantom>`) are ready, the experimental simulation can be readily done as following:

.. code-block:: lua
  
  local result = run{} 

The parameter structure is summarized as follow:

      .. list-table:: 
        :header-rows: 1
        :widths: 25 35 140

        * - Parameter
          - Mandatory/Optional
          - Content
        * - exp 
          - M
          - pulse sequence.
        * - spinsys
          - M
          - required by spectroscopy. 
        * - phantom 
          - M
          - required by imaging.
        * - supp
          - O
          - option for specifying user parameter list file for phantom.

.. code-block:: lua 
 
 result =run{exp =fid, spinsys =acrylic} -- spectroscopy
 result =run{exp =se, phantom ="circles.h5", supp ="config.txt"} -- imaging