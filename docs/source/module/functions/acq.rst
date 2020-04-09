:orphan:

***********
acq
***********
The data acquisition block ``acq`` was designed for signal observation, which simulate the quadrature detection of the receiver in MR scanner. E.g. 
     .. code-block:: lua 

       acq{np = 1024, sw = 1e4}

 performs 1H acquisition, in which ``sw`` specifies the spectral width in Hz, together with number of data points ``np``, the total acquisition time can be determined. 

     .. note::
        
        Similar to the r.f. blocks, explicit nuclear isotope should be appended to channel for other nuclei.


     .. code-block:: lua 

        acq{np = 1024, sw = 1e4, channel ="13C"}



It is worth emphasizing that Spin-Scenario is able to perform more complicated acquisition, **any of interested density states can be specified as the observer**. For example, 
     .. code-block:: lua 

       acq{np = 1024, sw = 1e4, observer = "2*I1zI2z"}
 can be used to acquire projection of the evolution state of spin system to the observation state 2I1zI2z.