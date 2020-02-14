:orphan:

*******
hardRF
*******

The most common used rectangle or hard pulses are realized by ``hardRF`` block.  
For example, 

     .. code-block:: lua 

       local rf1 = hardRF{beta =60, phase="-y"}

is a 60 degree -y-pulse for proton 1H excitation. The pulse width is proportional to that of the standard 90 degree pulse, which can be globally declared via ``pw90{}`` in advance. 

     .. code-block:: lua 

       pw90{5} -- 5 us for standard 90 degree pulse.

Alternatively, adding an explicit ``width`` key to the table will achieve a similar block. 

Besides, the valid phase options are not merely limited to the four quadrant axes, arbitrary angle in degree is also generally supported. 

For other nuclei, the explicit nuclear isotope 
such as ``channel = "13C"`` is required. Moreover, multi-channel excitation 
pulse can be achieved easily by extending both ``channel`` and ``phase`` respectively. 
For instance, 

     .. code-block:: lua 

       local rf2 = hardRF{beta =180, channel ="1H|13C", phase="x|x"}
       write("rf2.RF",rf2)

gives a refocusing pulse applied on a heteronuclear 1H-13C spin system.
