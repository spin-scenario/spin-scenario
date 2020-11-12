:orphan:

*******
hardRF
*******
The most common used rectangle or hard pulses are realized by **hardRF** block. 

Turn to :doc:`shapedRF<shapedRF>` for shaped pulses.

* **Create a pulse**
------------------------
  
  The syntax to create a rectangle pulse is quite simple:

  .. code-block:: lua 
        
        local rf = hardRF{}

  The parameter structure is summarized as follow:

  .. list-table:: 
    :header-rows: 1
    :widths: 25 35 140

    * - Parameter
      - Mandatory/Optional
      - Content
    * - beta
      - M
      - Flip angle in degree.
    * - channel
      - O
      - Nuclear isotope(s) for this pulse, default is ``"1H"``. For other nuclei, explicit nuclear isotope such as ``channel = "13C"`` or ``channel = "1H|13C"`` is required. Note for more channels, you only need to seperate the parameter for each channel with **|**.
    * - phase
      - O
      - Phase for this pulse, default is ``"x"`` phase. Shortcut symbols include ``"x"``, ``"-x"``, ``"y"``, ``"-y"``. You can also directly assign a numerical phase in degree such as ``"30"`` for customized phase. Similar to the `channel`, you only need to seperate the phase for each channel with **|**.
    * - width
      - O
      - Pulse duration in ms. If assigned, the amplitude will be determined automatically by the flip angle `beta`.


  .. note::
	  
    Pulse duration is proportional to that of the standard 90 degree pulse (default 5 microsecond), which indicates that the pulse amplitude is fixed for different flip angles. 
      
    The duration of standard 90 degree pulse can be globally declared as follow:

    .. code-block:: lua 

      pw90{10} -- unit in microsecond.

* **Demo script**
------------------------

  .. literalinclude:: ../../../../examples/seq/rf/hardRF.lua
	  :linenos:

