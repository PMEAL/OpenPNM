.. _bond_percolation_example:

===============================================================================
Standard Bond Percolation Simulations
===============================================================================
This examples illustrates how to set up an OpenPNM simulation to perform a typical bond and site percolation-type simulation.

As usual, start by importing OpenPNM, the necessary model libraries and retrieving the controller object:

.. code-block:: python

    import scipy as sp
    import OpenPNM
    from OpenPNM.Physics import models as pm
    from OpenPNM.Utrateilities.misc import tic, toc
    sim = OpenPNM.Base.Controller()

Next, generate a cubic Network of specified size that will contain the topological connections:

.. code-block:: python

    pn = OpenPNM.Network.Cubic(shape=[20,20,20],spacing=1)

In this case it is not necessary to create a Geometry object since we don't need pore and throat sizes.  Instead, create a Phase object that will act as our conductor:

.. code-block:: python

    ionomer = OpenPNM.Phases.GenericPhase(network=pn,name='ionomer')

This Phase object will provide the basic properties of the conducting phase, but we need to create a Physics object since the Algorithm object look there to find the necessary conductance values:

.. code-block:: python

    phys = OpenPNM.Physics.GenericPhysics(network=pn,phase=ionomer,pores=pn.Ps,throats=pn.Ts)
    phys['throat.random_seed'] = sp.rand(pn.Nt,)  # Assign random numbers to throats
    phys['pore.random_seed'] = sp.rand(pn.Np,)  # Assign random numbrers to pores
    phys['throat.ionic_conductance'] = 1.0

Finally, create an Algorithm object for performing the percolation calculations:

.. code-block:: python

    Seff = OpenPNM.Algorithms.OhmicConduction(network=pn,phase=ionomer)
    Seff.set_boundary_conditions(bctype='Dirichlet',bcvalue=1,pores=pn.pores('top'))
    Seff.set_boundary_conditions(bctype='Dirichlet',bcvalue=0,pores=pn.pores('bottom'))

For bond percolation, throats are progressively dropped from the Network and the overall Network conductivity is determined.  This is accomplished by running the Algorithm inside a loop, and to decreasing the number fraction of throat connections in the Network on each loop:

.. code-block:: python
    
    S = []  # Create an empty list to store the effective conductivity results
    phis = sp.logspace(-2,-0.01,30).tolist()  # Create a list of number fractions to simulate
    for phi in phis:
        Ts = phys['throat.random_seed']<phi  # Select
        phys['throat.ionic_conductance'][Ts] = 1.0e-5  # Set non-conducting throats to l
        Seff.run(conductance='throat.ionic_conductance',quantity='voltage')
        i = Seff.rate(outlets)
        S.append(float(i*L/((1-0)*A)))

The stored results can be visualized with Matplotlib:

.. code-block:: python

   plt.plot(phis[i],S[i])

To simulation site percolation, pores or sites are progressively removed as are all that pores neighboring throats.  Again this is done insdie a loop:

.. code-block:: python

    # Reset conductance values from bond percolation simulations
    phys['throat.conductance'][Ts] = 1.0
    S = []  # Create an empty list to store the effective conductivity results
    phis = sp.linspace(0.95,0.01,30).tolist()  # Create a list of number fractions to simulate
    for phi in phis:
        Ps = phys['pore.random_seed']>phi  # Select a fraction of pores in the network
        Ts = pn.find_neighbor_throats(pores=pn.toindices(Ps))  # Find the throats connected to pores
        phys['throat.conductance'][Ts] = 1.0e-5  # Set non-conducting throats to low value
        Seff.run(conductance='throat.conductance',quantity='voltage')
        i = Seff.rate(pn.pores('bottom'))
        S.append(i)
    

The stored results can be visualized with Matplotlib on the same axes as the bond percolation results:

.. code-block:: python

    plt.plot(phis,sp.log(S),'b.-')

The percolation threshold for bond percolation is lower that for site percolation.  The theoretical values are 24.88 and 31.16 `respectively <http://en.wikipedia.org/wiki/Percolation_threshold#Thresholds_on_3D_lattices>`_, which agrees with the present results considering the small size of the Network used here.  





