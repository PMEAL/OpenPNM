###############################################################################
CHANGE LOG
###############################################################################

===============================================================================
Version 1.2 (TBD)
===============================================================================
* Added ``find_nearest_pores`` to GenericNetwork for finding pspatially nearby pores regardless of whether they are topologically connected.



===============================================================================
Version 1.1 (April 9th, 2015)
===============================================================================
* Improved the interaction with models on each object.
* Introduced the Controller class, which provides high level oversight to all simulations.  The Controller allows saving and loading of simulations and entire workspaces.  
* Added the ability apply source and sink terms to the transport solvers.  The ``set_source_term`` method was added to GenericLinearTransport class and a variety of generic source term have been added to the available Physics models.


