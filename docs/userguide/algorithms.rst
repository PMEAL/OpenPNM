.. _algorithms:

===============================================================================
Algorithms
===============================================================================

.. warning::

    The wide variety of possible arguments and configurations for **Algorithms** makes it difficult to define a common behavior and requirements for initialization of this class.  In general, an **Algorithm** is instantiated by passing in a **Network** object, and the rest is up to the programmer of the class.  
	
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Basic Usage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
An **Algorithm** object can be instantiated by sending in a **Network** object:

>>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
>>> alg = OpenPNM.Algorithms.GenericAlgorithm(network=pn)

However, before the ``run`` method for this **Algorithm** can be called several steps need to be taken.  Most **Algorithms** have different requirements for what needs to be specified.  It is recommended to consult the docstring in each **Algorithm** for specific instructions.

-------------------------------------------------------------------------------
Adding Boundary Conditions
-------------------------------------------------------------------------------
The ``set_boundary_conditions`` method in the *GenericAlgorithm* class *can be* used to set boundary conditions on **Algorithm** objects.  Some algorithms are *not* required to use this method.  For instance OrdinaryPercolation algorithm simply requires a list of input pores, so a more complex approach is not needed (however it is possible).  In the case of transport simulations, however, this is quite essential:

>>> P1 = pn.pores('top')
>>> alg.set_boundary_conditions(bctype='Dirichlet',bcvalue=0.6,pores=P1, component=air)
>>> P2 = pn.pores('bottom')
>>> alg.set_boundary_conditions(bctype='Neumann',bcvalue=1e-7,pores=P2, component=air)

The ``set_boundary_conditions`` method does two things: it creates a *label* dictionary on the Algorithm object named according to the name of the component and also the type of boundary conditions applied, and applies that label to the specified locations. It also creates a *property* dictionary with the specified boundary values in the specified pores or throats.

For transport simulations, the *'bctype'* argument can be only one of the *reserved* keywords: *'Dirichlet'*, *'Neumann'* and *'Neumann_group'*. The first two keywords are self-explanatory, and the last one means that the total rate through *all* the specified locations has been set.

-------------------------------------------------------------------------------
Sharing Simulation Results
-------------------------------------------------------------------------------
Because each *Algorithm* object is Python *dictionary* it can and does internally store the results of any calculations it performs.  This prevents the objects in the simulation from becoming *polluted* with a large number of *property* arrays and prevent data from being overwritten.  In many situations, however, the results of one simulation are needed by another.  The *GenericAlgorithm* class provides a method called ``return_results`` which is not implemented; it is the responsibility of each **Algorithm** object to implement this method so that the call to ``alg.return_results(args)`` will send the important data to the correct places.  For instance, the **FickianDiffusion** algorithm will write *'pore.mole_fraction'* to the Phase object that it received as the *'phase'* argument.

An example of this can be found in the Function Reference for Ordinary Percolation.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Overview of Included Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
OpenPNM comes with several of the most common and widely used algorithms.  These include simulating transport phenomena in the pore space, ordinary percolation to simulate a drainage experiment, and invasion percolation to simulate fluid infiltration.

-------------------------------------------------------------------------------
Steady State Transport Algorithms
-------------------------------------------------------------------------------
By using pore scale physics, OpenPNM is capable of steady state simulation of Fickian diffusion, permeability calculation, heat and electron conduction through the phases in the network.  

In each of these algorithms, by using pipe network (or electrical resistor network) analogy, transport conductance of each conduit (a throat plus half of each adjoining pore) will be calculated.  Then, applying the conservation equation to each pore yields a sparse set of linear equations that can be solved with the appropriate boundary conditions to give the values of the desired quantity in each pore.


-------------------------------------------------------------------------------
FickianDiffusion
-------------------------------------------------------------------------------
It applies Fickian diffusion equation and uses the diffusive conductance of the desired phase to determine the mole fraction of that phase in each pore.  This binary diffusion algorithm, can be used for  unimolecular diffusion or equimolar counter diffusion.  In the case of unimolecular diffusion, however, the conversion of the mole fraction and BCs should take place outside of the algorithm.

-------------------------------------------------------------------------------
StokesFlow
-------------------------------------------------------------------------------
It applies a fluid flow transport equation such as Hagen-Poiseuille equation and uses the hydraulic conductance of the desired phase to determine the pressure of that phase in each pore.  

-------------------------------------------------------------------------------
OhmicConduction
-------------------------------------------------------------------------------
It applies Ohm's law to simulate electron or ion conduction and uses the electrical conductance of the desired phase to determine the voltage of that phase in each pore.  

-------------------------------------------------------------------------------
FourierConduction
-------------------------------------------------------------------------------
It applies Fourier equation to simulate heat conduction and uses the thermal conductance of the desired phase to determine the temperature of that phase in each pore.  

-------------------------------------------------------------------------------
OrdinaryPercolation
-------------------------------------------------------------------------------
OpenPNM includes a percolation algorithm for simulating drainage capillary pressure curves, where a specified pressure is applied and the volume of non-wetting liquid injected is monitored.  The most common version of this experiment is Mercury Intrusion Porosimetry, which is commonly used to measure pore size distributions.  The drainage process can be best simulated with a so-called 'access limited ordinary percolation' algorithm (ALOP).  ALOP is available as an additional option for this ordinary percolation algorithm.   

-------------------------------------------------------------------------------
InvasionPercolation
-------------------------------------------------------------------------------
Invasion percolation is related to OrdinaryPercolation, but it invades one pore at a time.  This algorithm is implemented in basic form in OpenPNM, using priority queues based on heaps.  This allows the calculation to reasonable fast for even large domains.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Creating Customized Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For description of how to create customized algorithms see :ref:`Customizing OpenPNM<customizing>`