:orphan:

***********
acq
***********
The data acquisition block **acq** was designed for signal observation, which simulate the quadrature detection of the receiver in MR scanner. 

   .. code-block:: lua 
        
        local adc = acq{}

   The parameter structure for **acq** is summarized as follow:

   .. list-table:: 
    :header-rows: 1
    :widths: 25 35 140

    * - Parameter
      - Mandatory/Optional
      - Content
    * - np
      - M
      - Number of data points.
    * - sw
      - M
      - Spectral width in Hz.
    * - channel
      - O
      - Nuclear isotope for this acquisition, default is ``"1H"``. For other nuclei, explicit nuclear isotope such as ``channel = "13C"`` is required. Currently only one channel acquisition is supported.
    * - observer
      - O
      - Spin-Scenario is able to perform more complicated acquisition, any of interested density states can be specified as the observer. E.g. ``observer = "2*I1zI2z"`` can be used to acquire projection of the evolution state of spin system to the observation state 2I1zI2z.

   .. note::
	  
      * Use **zero_padding{}** to add zeros to end of raw signal when needed.     
      .. code-block:: lua 
        
         zero_padding{4096}  -- Number of data points will be finally 4096.


      * Use **apodization{}** to perform signal apodization.     
      .. code-block:: lua 
        
         apodization{5}  -- Decay rate (default 10) for exp window function.