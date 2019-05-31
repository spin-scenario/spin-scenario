
*******
Utility
*******

Visualization
=============
 Spin-Scenario mainly use `Gnuplot <http://www.gnuplot.info/>`_ to visualize data.

Time–frequency analysis
=======================
This is an useful tool for the characteristics of shaped pulses.

signal from external files
--------------------------
.. code-block:: lua 

    specgram('wav.txt', {fs = 22050, window = 'hammning', wlen = 512, overlap = 0.75, nfft = 4096, style = 'dB'})

- ``fs``: sampling frequency in Hz.
- ``window``: window function such as ``'hammning'``, ``'gauss'``, etc. 
- ``wlen``: window length.
- ``overlap``: overlap ratio.
- ``nfft``: FFT number.
- ``style``: if assigned with ``'dB'``, the magnitude specgram will be shown in 20*log().

.. note::
	By default, the ``specgram`` function read only the 1st column data from the external file, which means the default signal is real signal. To analyze a complex signal (e.g. shaped pulses in ux/uy mode), you may need to specify the real and imag part with ``col = '1 2'``.
.. code-block:: lua

    specgram('coop_pulse.txt', {col ='6 7', fs = 219888, wlen = 128, overlap = 0.9, nfft = 2048})

signal within SSL
-----------------
For shaped pulses defined or optimized within the Spin-Scenario, you can analyze it directly like this:

.. code-block:: lua
    :linenos:

    local rf180 =shapedRF{pattern="sinc", width= 2, step=100}
    specgram(rf180, {wlen = 16, overlap = 0.9, nfft = 2048})
