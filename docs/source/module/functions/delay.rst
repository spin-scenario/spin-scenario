:orphan:

***********
delay
***********
The ``delay`` block is quite simple, there is only one parameter ``width`` (unit in **ms**) as below: 
     .. code-block:: lua 

       local tau = delay{width = 1.7857}

Note that it is also possible to generate multi-delay, which is commonly used MR spectroscopy. E.g.
     .. code-block:: lua 

       local tau = delay{ width = "1.7857:3:8"}
will generate 8 delays with equal intervals start from 1.7857 ms to 3 ms.
