.. _quick_start:

Getting Started
###############

It is relatively easy to get going with a quick simulation in OpenPNM.  In fact the following code block produces a mercury intrusion simulation in just 9 lines:

.. code-block::

   import openpnm as op
   Lc = 1e-4
   Nx, Ny, Nz = (10, 10, 10)
   pn = op.network.Cubic(shape=[Nx, Ny, Nz], spacing=Lc)
   geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps,  throats=pn.Ts)
   Hg = op.phases.Mercury(network=pn)
   phys = op.physics.Standard(network=pn, phase=Hg, geometry=geo)
   mip = op.algorithms.Porosimetry(network=pn)
   mip.setup(phase=Hg)
   mip.set_inlets(pores=pn.pores(['left', 'right', 'top', 'bottom']))
   mip.run()

The results can be visualized with The ``mip.plot_intrusion_curve()`` giving something like this:

.. image:: https://user-images.githubusercontent.com/14086031/77930201-96363b80-7278-11ea-95fd-4a55fb1d6148.png
   :width: 800px

As another example, the permeability coefficient can be found as follows:

.. code-block::

   water = op.phases.Water(network=pn)
   phys = op.physics.Standard(network=pn, phase=water, geometry=geo)
   sf = op.algorithms.StokesFlow(network=pn, phase=water)
   Pin, Pout = (200000, 101325)
   sf.set_value_BC(pores=pn.pores('left'), values=Pin)
   sf.set_value_BC(pores=pn.pores('right'), values=Pout)
   sf.run()

The total flow rate into the domain through the boundary pores can be found using ``sf.rate(pores=pn.pores('left'))``.  The permeability coefficient can be found by inserting known values into Darcy's law as follows:

.. code-block:: python

   Q = sf.rate(pores=pn.pores('left'))
   A = Ny*Nz*Lc**2
   L = Nx*Lc
   mu = water['pore.viscosity'].mean()
   K = Q*mu*L/(A*(Pin-Pout))

It's also worth explaining how to adjust the pore size distribution of the network, so that the capillary curve and permeability coefficient can be changed to match known values.  The ``geo`` object controls the geometric properties, and it possess models to calculate values on demand. Let's change the pore size distribution to a Weibull distribution, but first let's store the existing values in a dummy variable so we can compare later.

.. code-block:: python

   geo['pore.old_diameter'] = geo.pop('pore.diameter')
   geo.add_model(propname='pore.diameter',
                 model=op.models.geometry.pore_size.weibull,
                 shape=0.5, loc=0, scale=1e-5)


Now you can run ``geo.show_hist(['pore.old_diameter', 'pore.diameter'])`` to get a quick glance at the histograms of the two distributions.

More complex tasks are explained in the `online examples <https://github.com/PMEAL/OpenPNM/tree/dev/examples>`_.
