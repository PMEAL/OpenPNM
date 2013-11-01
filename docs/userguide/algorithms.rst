. _algorithms:

###############################################################################
Algorithms
###############################################################################

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Percolation Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

===============================================================================
Invasion Percolation
===============================================================================
OpenPNM employs a sophisticated drainage algorithm that allows for a wide variety of initial conditions. The algorithm follows the basic tenants of the invasion percolation scheme, which was first described here [1]_, with a few clarifications.

1) Multiple, independent clusters of the invading phase can be simulated, originating from any pore, or group of pores in the network.
2) Clusters grow simultaneously. Pore invasion times are calculated from cluster volumes and prescribed cluster growth rates.
3) A list of pore indices must be specified as outlets. A simulation only reaches its breakthrough condition once each invading cluster has either reached one of these pores or has coalesced with a cluster that has.
4) Isolated, or trapped, defending phase is not considered. In a future release, trapping logic will be available.

*Initial Setup*
Two parameters must be specified: inlets and outlets. inlets is a list of arrays, where the arrays should be non-overlapping sets of available pore names. Each set is assumed to be a continuous cluster of invading fluid.
By default, the program assumes that clusters initially grow at initial rates, but the user may specify individual cluster growth rates. The clusters each begin to grow throughout the network once enough simulation time has passed to fill the original pore volumes.

*Cluster Volume Calculation*
At each pore filling step, the subsequent filling step is anticipated by searching the cluster's interface for the lowest capillary pressure barrier. It is assumed that the entire cluster will reach the pressure of this capillary barrier before the cluster can advance (quasi-static flow). To reach this pressure, not only should all invaded pores be completely filled, but each interfacial meniscus should inflate to the radius of curvature associated with the capillary pressure.  Ideally, the shape of the throat, the local fluid/fluid surface tension, and a local contact angle should all factor in to the calculation of each inflated menisci' contribution to the cluster volume. However, for this version of the algorithm, the approximation is made that meniscus volume is a linear function of cluster capillary pressure, intersecting zero volume at zero pressure and the volume of a hemisphere with the throat's diameter at the throat's barrier capillary pressure.

.. [1] D Wilkinson and J F Willemsen 1983 J. Phys. A: Math. Gen. 16 3365 http://dx.doi.org/10.1088/0305-4470/16/14/028

===============================================================================
Access Limited Ordinary Percolation
===============================================================================
OpenPNM also includes an algorithm for simulating drainage capillary pressure curves, where a specified pressure is applied and the volume of non-wetting liquid injected is monitored.  The most common version of this experiment is Mercury Intrusion Porosimetry, which is commonly used to measure pore size distributions.  The drainage process is best simulated with a so-called 'access limited ordinary percolation' algorithm (ALOP).  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Transport Algorithms
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
OpenPNM is capable of simulating mass transfer diffusion, permeability calculation, heat and electron conduction through the fluids in the network, by using pore scale physics which is described here [2]_, .

To execute a transport algorithm for a fluid on a network (each with known recipe), it is recommended to: 

1- Create an algorithm object.

2- Apply boundary conditions to the object.

3- Run the algorithm object for the active fluid in the network.

*Boundary Conditions:*
For boundary conditions, algorithm needs to know values of boundary conditions, their types (Dirichlet, Neumann, etc)  and also the pore numbers. Therefore, for any network applying boundary conditions to arbitrary pores (whether they are internal or boundary pores) includes two steps:

1- Creating two 1D zero arrays of length Np, which Np is the total number of the pores (internal and boundary pores)

2- Assigning boundary types and boundary values to the pores as BCtypes and BCvalues, respectively.

The followings are useful hints for applying boundary conditions:

- BCtypes and BCvalues are bound to the object algorithm.
- Types of boundary condition: Dirichlet = 1, Neumann_flux = 2, Neumann_insulated = 3, Neumann_rate = 4
- In Fickian algorithm, positive value for *Neumann_rate* or *Neumann_flux* for a pore means that the quantity of interest leaves the pore, but for any other algorithms, positive Neumann value for a pore means that the quantity of interest enters this pore. This is because of variable transformation in Fickian algorithm from mole fraction of active fluid to natural logarithm of stagnant film mole fraction.
- For a pore with *Neumann_rate* type, the boundary value of the pore might not represent the rate of quantity for just that pore, but can imply that this pore belongs to a group of the pores which this amount of rate enters/leaves all of them. By this assumption, when the rate for each individual pore is unknown, it is still possible to use *Neumann_rate* boundary conditions. However the user can always apply the rate for just a single pore, but this rate should be unique among all of the pores in the network. In other words, if two pores with *Neumann_rate* boundary type, have different values, algorithm object considers them as their individual rate, but if their rate is exactly the same, algorithm will assume that this rate is the total rate which enters/leaves both of them.  
- *Neumann_insulated* is equivalent to Neumann_flux boundary condition when flux is zero. Therefore, there is no need to define BCvalues for this kind of boundary condition.
- For a boundary pore in the network (a pore in boundary faces), if the type of boundary condition has not been specified, it is assumed to be *Neumann_insulated* type. This will help user to apply boundary conditions with less lines of code.
- Units for all values should be in SI.
- Each location in BCtypes and BCvalues arrays will determine the pore number.

For example, if BCtypes and BCvalues are defined as:

	BCtypes = [2, 1, 4 ,0, 4, 3, 4, 1]

	BCvalues = [0.1, 0.5, 0.87, 0, -0.35, 0, 0.87, 0.30]
 
It means that:

for pore 0: Neumann, flux = 0.1

for pore 1: Dirichlet, value = 0.5

for pore 2: Neumann, rate = 0.0087 (hint: Since there are two pores (2,6) with *Neumann_rate* type which have the exact same amount of rate, algorithm assumes that 0.0087 is the rate of quantity of interest which enters/leaves both pore 2 and 6)

for pore 3: Internal pore without imposed boundary condition (hint: If pore 3 is a boundary pore (a pore in a boundary face), algorithm by default assumes that, this is *Neumann_insulated* pore.)

for pore 4: Neumann, rate= -0.35 (hint: There is only one pore with Neumann_rate type and value of -0.35. So, algorithm assumes that 0.35 is the rate of quantity of interest which is only entering/leaving pore 4)

for pore 5: Neumann_insulated, value=0

for pore 6 : Neumann, rate=0.0087 (hint: Refer to pore 2)

for pore 7 : Dirichlet, value = 0.30

*Example for Applying boundary conditions in OpenPNM:*
For the network ``pn``, applying Dirichlet boundary condition to face 1 and face 6 with values of 8e-2 and 8e-1, can be achieved as follows:

.. code-block:: python

   >>> BCtypes = sp.zeros(pn.get_num_pores())
   >>> BCvalues = sp.zeros(pn.get_num_pores())
   >>> BCtypes[pn.pore_properties['type']==1] = 1
   >>> BCtypes[pn.pore_properties['type']==6] = 1
   >>> BCvalues[pn.pore_properties['type']==1] = 8e-2
   >>> BCvalues[pn.pore_properties['type']==6] = 8e-1

As can be seen, by this method applying boundary conditions for any type of the network and to any pore in the network is straightforward.


 
.. [2]  Gostick, J.T., et al., Pore network modeling of fibrous gas diffusion layers for polymer electrolyte membrane fuel cells. Journal of Power Sources, 2007. 173(1): p. 277-290.


===============================================================================
Diffusion
===============================================================================

-------------------------------------------------------------------------------
Fickian
-------------------------------------------------------------------------------
OpenPNM applies Fickian diffusion equation to determine flux between the pores. By using pipe network (or electrical resistor network) analogy, conductance of each conduit (a throat plus half of each adjoining pore) will be calculated. Diffusive conductance for a fluid in each conduit depends on the geometry of the conduit, the binary diffusion coefficient, the total molar concentration and also occupancy of the fluids in the bond.

Using Fickian equation for diffusion of fluid ``A`` through stagnant fluid ``B``, and then applying species conservation equation to each pore in the network yields a sparse set of linear equations that can be solved with the appropriate boundary conditions to give mole fraction of the active fluid in each pore.

To run Fickian algorithm for ``active_fluid = air`` in the network ``pn`` :

.. code-block:: python

   >>> Fickian_alg = OpenPNM.Algorithms.FickianDiffusion()
   >>> Fickian_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
   >>> Fickian_alg.run(pn,active_fluid=air)


This algorithm will store mole fractions of ``air`` in ``air.pore_conditions['mole_fraction']``. Conduit conductance can also be found in ``air.throat_conditions['diffusive_conductance']``. Effective diffusivity between face ``i`` and face ``j`` in a cubic network, can be achieved as follows:


.. code-block:: python

   >>> Fickian_alg.calc_eff_diffusivity_cubic(face1=i,face2=j)

   
-------------------------------------------------------------------------------
Stefan-Maxwell (TODO)
-------------------------------------------------------------------------------

===============================================================================
Permeability
===============================================================================

-------------------------------------------------------------------------------
Hagen-Poiseuille
-------------------------------------------------------------------------------

The fluid flow transport  in the network based on Hagen-Poiseuille equation can be found in the same manner as Fickian diffusion algorithm. For running permeability algorithm, we can write:

.. code-block:: python

   >>> permeability_alg = OpenPNM.Algorithms.Permeability()
   >>> permeability_alg.set_boundary_conditions(types=BCtypes,values=BCvalues)
   >>> permeability_alg.run(pn,active_fluid=air)

Permeability algorithm stores pressure gradient in the network for active fluid in  ``air.pore_conditions['pressure']`` and conduit conductance in ``air.throat_conditions['hydraulic_conductance']``. 
Effective permeability between face ``i`` and face ``j`` in a cubic network, can be achieved as follows:

.. code-block:: python

   >>> permeability_alg.calc_eff_permeability_cubic(face1=i,face2=j)



-------------------------------------------------------------------------------
Klinkenburger Slip-flow (TODO)
-------------------------------------------------------------------------------


===============================================================================
Filtration (TODO)
===============================================================================
TODO

===============================================================================
Bubble Growth and Condensation (TODO)
===============================================================================

===============================================================================
Writing Custom Algorithms
===============================================================================
- See developers guide
