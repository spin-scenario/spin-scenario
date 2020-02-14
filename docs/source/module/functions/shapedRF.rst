:orphan:

*********
shapedRF
*********
Shaped pulses has been widely used to increase excitation bandwidth, achieve desired profile over the entire bandwidth, and improve polarization and coherence transfer efficiencies in numerous NMR experiments. To this end, ``shapedRF`` was provided as a general interface for routine patterns (such as Sinc, Gaussian, Rectangle, etc.), external shape files and complex mathematical expressions. 

  .. literalinclude:: ../../../../examples/seq/rf/shapedRF.lua
	  :linenos:

  |rf_sinc|

  |rf_shape_file|

  |rf_expr|

.. |rf_sinc| image:: ../../media/seq/rf_sinc.png
  :height: 480
  :align: middle

.. |rf_shape_file| image:: ../../media/seq/rf_shape_file.png
  :height: 480
  :align: middle

.. |rf_expr| image:: ../../media/seq/rf_expr.png
  :height: 480
  :align: middle  