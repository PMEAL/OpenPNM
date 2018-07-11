.. _percolation_guide:

================================================================================
Percolation
================================================================================

.. contents:: Page Contents
    :depth: 3

--------------------------------------------------------------------------------
Introduction
--------------------------------------------------------------------------------
In general terms percolation is transport associated with some threshold behaviour. Typically we are concerned with modelling
multiphase flow in porous media dominated by capillary forces and the threshold is the capillary pressure required
for a phase to enter a pore or throat by displacing another phase. However, other physical processes can be considered Percolation
such as the ising model. Percolation is associated with the presence of phases under certain conditions and the connectivity of these phases.
When referring to transport of wetting and non-wetting fluids in porous media, drainage is defined as the displacement of the wetting phase by invasion of the non-wetting phase.
A non-wetting phase requires positive pressure to overcome surface tension and push out the wetting phase and so as the imposed pressure in the non-wetting phase increases,
it pushes further into the porous media. Capillary pressure is determined by the wettability of the phases w.r.t the solid material characterized by the contact angle, and the geometry of the material which together shape
the meniscus. A higher curvature requires more pressure and so a non-wetting phase will require greater pressure to squeeze into smaller spaces.

--------------------------------------------------------------------------------
Ordinary Percolation
--------------------------------------------------------------------------------

Ordinary Percolation (OP) has two modes: site and bond. In OpenPNM we refer to sites as pores and bonds as throats and we can refer to both generally as elements.
OP identifies clusters of elements that are connected and occupied by the same phase given the threshold pressure and the imposed pressure.
OpenPNM uses physical models to determine what the highest pressure is required to enter a pore or throat known as the entry pressure.
So to generate a porosimetry curve which is a sum of the volume of pore space occupied by the phases
of interest for a given pressure, we define a range of pressures and for each one we compare the value to the entry pressures of the elements (pores or throats, not both as discussed later).
If the entry pressure is greater than the current threshold value this element is considered invaded and can form part of an invading cluster. However, another step is required to determine whether
each cluster has access to an inlet. Invasion of the non-wetting phase must progress from an inlet and this is referred to as access limited. Another case is invasion of the wetting phase, referred to
as imbibition and this may progress by different physical mechanisms such as film growth which may be access unlimited as wetting films can permeate the entire network.
OP is a quasi-static algorithm which has more of a basis in graph theory than real physical simulations of transport in porous media. It is useful for gathering information about the pore size distribution,
but not really for simulating multiphysics. However, it's advantages is that it can be very fast and appropriate for simulating common experimental data such as mercury intrusion porosimetry.

--------------------------------------------------------------------------------
Invasion Percolation
--------------------------------------------------------------------------------

Invasion Percolation (IP) again typically has two modes: site and bond, but also allows for site and bond percolation which we describe in the next section.
Similarly to OP we are again concerned with analysis of the entry pressure. However, instead of identifying connected clusters and invading them all in one step, we identify the neighboring elements
of the invading cluster and invaded one at a time based on the path of least resistance. This method allows for a more accurate representation of transient flow and for more physical models associated with
the position and advancement of the meniscus within a given element. Phenomena such as trapping where clusters can become isolated, co-operative pore filling and snap off are also only possible with IP.
It is possible to define multiple inlet clusters which may progress at different rates and pressures, again allowing for more physical situations to be simulated. The draw-back to IP is that for larger networks
it can be significantly slower, although care has been taken to optimize the algorithms as much as possible. When running IP in site mode the capillary pressure of the pores are used and all throats connected to an
invaded pore are also considered to be invaded. Conversely, when running in bond mode, the entry pressure of the throats is used and connected pores are automatically invaded.
This is really a convention used to speed up calculations with reasoning being that throats are typically smaller than pores. Therefore, for drainage the throats require a higher capillary pressure and so once the meniscus has reached this point
it can freely enter a larger connected space making bond percolation the most appropriate. The reverse scenario is imbibition where larger spaces provide greater resistance to flow (the ink bottle effect) and so site percolation is appropriate.

--------------------------------------------------------------------------------
Mixed Invasion Percolation
--------------------------------------------------------------------------------

As mentioned before, Mixed Invasion Percolation, is a special case of IP where both pores and throats can be invaded on an individual basis, this is appropriate when the wettability of the invading and defending phases are similar,
in this case the porous media is said to have neutral wettability. Other factors other than simple pore and throat sizes can determine the shape and displacement of the meniscus and Mixed IP allows for processes in both pores and throats to happen in the same simulation such
as cooperative pore filling and throat snap-off.
