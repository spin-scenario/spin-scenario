
***********
Spin System
***********

The formation of spin systems in *Spin-Scenario* is based on quantum mechanics instead of on the semiclassical Bloch approach. Basically, the quantum state of a spin system can be represented by the density matrix or the density operator, which is essentially built based on irreducible spherical tensors. The theory details of the density operator are well described `here <https://www.sciencedirect.com/science/article/pii/S1090780705001564>`_.

spin system creation
------------------------
In  *Spin-Scenario*, a spin system can be easily created with a Lua table, including all necessary parameters such as ``spin``, ``zeeman``, ``jcoupling`` and ``relaxation``. For example, the 13C labelled alanine can be formed with following code snippet:

  .. code-block:: lua 

        local sys = spin_system{
        spin = "13C 13C 13C",
        zeeman = "1 scalar 15.74 kHz 3 scalar  -4.3 kHz", 
        jcoupling ="1 2 scalar 54.2 Hz 2 3 scalar 35.4 Hz"
        }

operator and state
---------------------
In general, the density operators or state operators can be represented as products of individual spin operators. The ``op`` and ``state`` functions are used for representing the density operator and the state operator respectively. The validated strings for each spin are 

    * I+
    * I-
    * Iz
    * Ie
E.g. a state I1pI3m can be defined as
  .. code-block:: lua 
        
        sys:state{1, "I+", 3, "I-"}

Moreover, taking advantage of the computer algebra system `YACAS <http://www.yacas.org/>`_,  Spin-Scenario also provides a very practical feature for the construction of more complex operators by ``expr_state`` that are a linear combination of arbitrary either spherical or cartesian operators. E.g. 
  .. code-block:: lua 

    sys:expr_state("2*I1xI3x-0.5*I2xI3y")

is a straightforward generation of the combined state 2I1xI3x-0.5I2xI3y. The Yacas integration was fully utilized throughout Spin-Scenario to offer a nature manner for kinds of expression evaluation.

For a specific spin system, it is also possible to retrieve its free Hamiltonian, total Hamiltonian and the relaxation superoperator via ``sys.L0``, ``sys.L`` and ``sys.R`` respectively.
