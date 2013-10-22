###############################################################################
Fluids and Phases Module
###############################################################################
Virtually all algorithms in OpenPNM require fluid properties of some sort, either directly or through the pore scale physics methods.  The Fluids module in OpenPNM provide a simple way to ensure that all fluid properties are available to the network when needed, and that they are up to date with changes in the network conditions.  

===============================================================================
What is a Fluid in OpenPNM?
===============================================================================
A fluid object in OpenPNM does **not** contain the properties of the fluid.  Instead it contains the *recipe* for generating the fluid's properties.  The reason for this rather unconventional definition arises from the following fact:  in order for a fluid to have properties, it must have (at least) a temperature and a pressure.  Any temperature and pressure assigned to the fluid upon its creation (in order to define properties) would be totally arbitrary.  Instead, OpenPNM leaves the fluid in an undefined state until it is *assigned* a pore network object, which **does** have a temperature and pressure.  Once the association is made the fluid properties can be determined by executing it's *recipe*.  

This approach is not just philosophical.  The process of *assigning* a fluid to a pore network in stores the fluid *recipe* on the pore network under `pn.phases['fluid_name']`.  This means that the pore network can then update or recalculate the fluid properties whenever an algorithm changes some environmental condition of the network.  The details of these steps are described below.  

-------------------------------------------------------------------------------
Creating a Fluid
-------------------------------------------------------------------------------



-------------------------------------------------------------------------------
Assigning a Fluid to a Network
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------
Refreshing a Fluid for Changes in Network Conditions
-------------------------------------------------------------------------------