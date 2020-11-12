
***************
Optimal Control
***************

The primary goal of optimal control of spin system is usually to maximize fidelity to desired state or unitary using optimized pulse shape, with additional constraints associated with specific experimental requirements (e.g., B0 and rf inhomogeneity, rf power limitation). There have been numerous methods for this purpose, and the algorithms involved are mostly based on gradient methods, such as `GRAPE <http://refhub.elsevier.com/S1090-7807(19)30022-9/h0105>`_ and `Krotov <http://refhub.elsevier.com/S1090-7807(19)30022-9/h0110>`_. 

Spin-Scenario offers the following two optimal control schemes which enables an efficient optimization of shaped pulses or `multiple cooperative pulses <http://refhub.elsevier.com/S1090-7807(19)30022-9/h0145>`_.
 
* :doc:`rf_optimizer<functions/rf_optimizer>`, for single shaped pulse optimization.
* :doc:`multi_rf_optimizer<functions/multi_rf_optimizer>`, for cooperative pulses optimization.
