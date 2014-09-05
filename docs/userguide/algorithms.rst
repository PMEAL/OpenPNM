.. _algorithms:

===============================================================================
Algorithms
===============================================================================


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

.. warning::

    The wide variety of possible arguments and configurations for algorithms makes it difficult to define a common behavior and requirements for this initialization of this class.  In general, an Algorithm is instantiated by passing in a Network object, and the rest is up to the programmer of the class.  

An Algorithm object can be instantiated by sending in a Network object:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
>>> alg = OpenPNM.Algorithms.GenericAlgorithm(network=pn)

However, before the ``run`` method for this Algorithm can be called several steps need to be taken.  Most Algorithms have different requirements for what needs to be specified.  It is recommended to consult the docstring in each Algorithm for specific instructions.

-------------------------------------------------------------------------------
Adding Boundary Conditions
-------------------------------------------------------------------------------
The ``set_boundary_conditions`` method in the GenericAlgorithm class *can be* used to set boundary conditions on Algorithms.  Algorithms are *not* required to use this method.  For instance OrdinaryPercolation algorithm simply requires a list of input pores so a more complex approach is not needed.  In the case of transport simulations, however, this is quite essential:

>>> P1 = pn.pores('top')
>>> alg.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=P1)
>>> P2 = pn.pores('bottom')
>>> alg.set_boundary_conditions(bctype='Neumann',bcvalue=1e-7,pores=P2)

The ``set_boundary_conditions`` method does 2 things: it creates a 'label' dictionary on the Algorithm object named according the type of boundary conditions applied, and applies that label to the specified pores, and it creates a 'property' dictionary with the specified boundary values in the specified.

The 'bctype' argument accepts 'reserved' keywords: 'Dirichlet', 'Neumann' and 'Neumann_group'.  The fist two should be self-explanatory, and the last one means that the total flux through 'all' the specified pores set.

-------------------------------------------------------------------------------
Sharing Simulation Results
-------------------------------------------------------------------------------
Because each Algorithm object is a ``dict`` it can internally store the results of any calculations it performs.  This prevents the objects in the simulation from becoming 'polluted' with a large number of 'property' arrays and prevent data from being overwritten.  In many situations, however, the results of one simulation are needed by another.  The GenericAlgorithm class provides a method called ``update_results`` which is not implemented.  It is the responsibility of each Algorithm object to implement this method so that the call to ``alg.update(args)`` will send the important data to the correct places.  For instance, the FickianDiffusion algorithm will write 'pore.mole_fraction' to the Phase object that it received as the 'active_phase' argument.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Overview of Included Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
OpenPNM comes with several of the most common and widely used algorithms.  These include simulating transport phenomena in the pores space, ordinary percolation to simulate a drainage experiment, and invasion percolation to simulate fluid infiltration.

-------------------------------------------------------------------------------
FickianDiffusion
-------------------------------------------------------------------------------
in progress

-------------------------------------------------------------------------------
StokesFlow
-------------------------------------------------------------------------------
in progress

-------------------------------------------------------------------------------
OhmicConduction
-------------------------------------------------------------------------------
in progress

-------------------------------------------------------------------------------
FourierConduction
-------------------------------------------------------------------------------
in progress

-------------------------------------------------------------------------------
OrdinaryPercolation
-------------------------------------------------------------------------------
in progress

-------------------------------------------------------------------------------
InvasionPercolation
-------------------------------------------------------------------------------
in progress

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Creating Customized Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized algorithms see :ref:`Customizing OpenPNM<customizing>`