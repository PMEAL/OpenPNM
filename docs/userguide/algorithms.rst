. _algorithms:

###############################################################################
Algorithms
###############################################################################

===============================================================================
Percolation Algorithms
===============================================================================

-------------------------------------------------------------------------------
Invasion Percolation
-------------------------------------------------------------------------------
OpenPNM employs a sophisticated drainage algorithm that allows for a wide 
variety of initial conditions. The algorithm follows the basic tenants of the 
invasion percolation scheme, which was first described here [1]_, with a few 
clarifications.

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

-------------------------------------------------------------------------------
Access Limited Ordinary Percolation
-------------------------------------------------------------------------------
OpenPNM also includes an algorithm for simulating drainage capillary pressure curves, where a specified pressure is applied and the volume of non-wetting liquid injected is monitored.  The most common version of this experiment is Mercury Intrusion Porosimetry, which is commonly used to measure pore size distributions.  The drainage process is best simulated with a so-called 'access limited ordinary percolation' algorithm (ALOP).  

===============================================================================
Transport Algorithms===============================================================================This documentation is being rewritten, sorry for the inconvenience.


+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Diffusion
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

-------------------------------------------------------------------------------
FickianDiffusion
-------------------------------------------------------------------------------
OpenPNM applies Fickian diffusion equation to determine flux between the pores.

-------------------------------------------------------------------------------
StefanMaxwell (TODO)
-------------------------------------------------------------------------------
This documentation is being rewritten, sorry for the inconvenience.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Permeability
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
This documentation is being rewritten, sorry for the inconvenience.

-------------------------------------------------------------------------------
StokesFlow
-------------------------------------------------------------------------------
This documentation is being rewritten, sorry for the inconvenience.

-------------------------------------------------------------------------------
Klinkenburger Slip-flow (TODO)
-------------------------------------------------------------------------------
This documentation is being rewritten, sorry for the inconvenience.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Filtration (TODO)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TODO

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Bubble Growth and Condensation (TODO)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TODO

===============================================================================
Writing Custom Algorithms
===============================================================================
TODO

