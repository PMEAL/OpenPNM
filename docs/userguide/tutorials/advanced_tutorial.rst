.. _advanced_tutorial3

###############################################################################
Tutorial 3 of 3: Advanced Topics and Usage
###############################################################################

.. contents:: Topics Covered in this Tutorial

**Learning Outcomes**

#. Use different methods to add boundary pores to a network
#. Manipulate network topology by adding and removing pores and throats
#. Explore the ModelsDict design, including copying models between objects, and changing model parameters
#. Write a custom pore-scale model and a custom Phase
#. Access and manipulate objects associated with the network
#. Combine multiple algorithms to predict relative permeability
#. Use the workspace manager to save and load, clone and purge simulations

===============================================================================
Build and Manipulate Network Topology
===============================================================================

For the present tutorial, we'll keep the topology simple to help keep the focus on other aspects of OpenPNM.

.. code-block:: python

    >>> import scipy as sp
    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 10], spacing=0.00006, name='net')

-------------------------------------------------------------------------------
Adding Boundary Pores
-------------------------------------------------------------------------------

When performing transport simulations it is often useful to have 'boundary' pores attached to the surface(s) of the network where boundary conditions can be applied.  The **Cubic** class has two methods for doing this: ``add_boundaries`` and ``add_boundary_pores``.  The first method automatically adds boundary to ALL six faces of the network and offsets them from the network by 1/2 of the value provided as the network ``spacing``.  The second method provides total control over which boundary pores are created and where they are positioned, but it more cumbersome to use.  Let's explore these two options:

.. code-block:: python

    >>> pn.Np
    1000
    >>> pn.Nt
    2700
    >>> pn.add_boundaries()
    >>> pn.Np
    1600
    >>> pn.Nt
    3300

Let's remove all these newly created boundary pores.  When they are created these pores are all automatically labeled with a label such as ``'top_boundary'``, so we can select all boundary pores by using the 'wildcard' feature in the ``pores`` look-up method to find all pores with a label containing the word ``boundary``.

.. code-block:: python

    >>> Ps = pn.pores('*boundary')  # Using the * wildcard

We can then ``trim`` these pores from the network using:

.. code-block:: python

    >>> pn.trim(pores=Ps)
    >>> pn.Np
    1000
    >>> pn.Nt
    2700

Note that all throats connecting to the trimmed pores were automatically removed since OpenPNM does not allow 'dangling' or 'headless' throats.

Now that ``pn`` is back to its original size, let's explore the second approach to apply boundary pores.

.. code-block:: python

    >>> Ps = pn.pores('top')  # Select pores on top of network
    >>> pn.add_boundary_pores(pores=Ps, offset=[0, 0, 0.00003],
    ...                       apply_label='top_boundary')
    >>> Ps = pn.pores('bottom')  # Select pores on bottom of network
    >>> pn.add_boundary_pores(pores=Ps, offset=[0, 0, -0.00003],
    ...                       apply_label='bottom_boundary')
    >>> pn.Np
    1200
    >>> pn.Nt
    2900

This approach requires more typing than the ``add_boundaries`` method, but allows for much finer control over how boundaries are created.

-------------------------------------------------------------------------------
Adding and Removing Pores and Throats
-------------------------------------------------------------------------------

OpenPNM uses a `list-based data storage scheme <topology>`_ for all properties, including topological connections.  One of the benefits of this approach is that adding and removing pores and throats from the network is essentially as simple as adding or removing rows from the data arrays.  The one exception to this 'simplicity' is that the ``'throat.conns'`` array must be treated carefully when trimming pores, so OpenPNM provides the ``extend`` and ``trim`` functions for added and removing, respectively.  To demonstrate, let's reduce the coordination number of the network to create a more random structure:

.. code-block:: python

    >>> Ts = sp.rand(pn.Nt) < 0.1  # Create a mask with ~10% of throats labeled True
    >>> pn.trim(throats=Ts)  # Use mask to indicate which throats to trim

When the ``trim`` function is called, it automatically checks the health of the network afterwards, so logger messages might appear on the command line if problems were found such as isolated clusters of pores or pores with no throats.  This health check is performed by calling the **Network**'s' ``check_network_health`` method which returns a **HealthDict** containing the results of the checks:

.. code-block:: python

    >>> a = pn.check_network_health()
    >>> pn.trim(pores=a['trim_pores'])

The **HealthDict** contains several lists including things like duplicate throats and isolated pores, but also a suggestion of which pores to trim to return the network to a healthy state.  Also, the **HealthDict** has a ``health`` attribute that is ``False``` is any checks fail.

===============================================================================
Define Geometry Objects
===============================================================================

The boundary pores we've added to the network should be treated a little bit differently.  Specifically, they should have no volume or length (as they are not physically representative of real pores).  To do this, we create two separate **Geometry** objects, one for internal pores and one for the boundaries:

.. code-block:: python

    >>> Ps = pn.pores('*boundary', mode='not')
    >>> geom = op.Geometry.Stick_and_Ball(network=pn, pores=Ps, throats=pn.Ts,
    ...                                   name='internal')
    >>> Ps = pn.pores('*boundary')
    >>> boun = op.Geometry.GenericGeometry(network=pn, pores=Ps, name='boundary')

The **Stick_and_Ball** class is preloaded with the pore-scale models to calculate all the necessary size information (pore diameter, throat lengths, etc).  The **GenericGeometry** class used for the boundary pores is empty and requires work:

.. code-block:: python

    >>> boun['pore.diameter'] = 0
    >>> boun['pore.volume'] = 0

These models are required for the Hagan-Poiseuille model. Most of them are straight-forward geometry calculations, except for the model used for ``'throat.diameter'``.  In this case the model looks into the neighbor pores, retrieves the two ``'pore.diameter'`` and uses the ``'max'`` value.  Because we set the boundary pores to have 0 diameter, this will naturally find result in the throat being assigned the diameter of the internal pore.

===============================================================================
Define Multiple Phase Objects
===============================================================================

In order to simulate relative permeability of air through a partially water-filled network, we need to create each **Phase** object.  OpenPNM includes pre-defined classes for each of these common fluids:

.. code-block:: python

    >>> air = op.Phases.Air(network=pn)
    >>> water = op.Phases.Water(network=pn)
    >>> water['throat.contact_angle'] = 110
    >>> water['throat.surface_tension'] = 0.072

-------------------------------------------------------------------------------
Aside: Creating a Custom Phase Class
-------------------------------------------------------------------------------

In many cases you will want to create your own fluid, such as an oil or brine, which may be commonly used in your research.  OpenPNM cannot predict all the possible scenarios, but luckily it is easy to create a custom **Phase** class as follows:

.. code-block:: Python
    :linenos:
    :caption: **Example of a Subclassed Phase**

    from OpenPNM.Phases import GenericPhase, models

    class Oil(GenericPhase):
        def __init__(self, **kwargs):
            super().__init__(**kwargs)
            self.models.add(propname='pore.viscosity',
                            model=models.misc.polynomial,
                            poreprop='pore.temperature',
                            a=[1.82082e-2, 6.51E-04, -3.48E-7, 1.11E-10])
            self['pore.molecular_weight'] = 116  # g/mol

* Creating a **Phase** class basically involves placing a series of ``self.models.add`` commands within the ``__init__`` section of the class definition.  This means that when the class is instantiated, all the models are added to *itself* (i.e. ``self``).

* ``**kwargs`` is a Python trick that captures all arguments in a *dict* called ``kwargs`` and passes them to another function that may need them.  In this case they are passed to the ``__init__`` method of **Oil**'s parent by the ``super`` function.  Specifically, things like ``name`` and ``network`` are expected.

* The above code block also stores the molecular weight of the oil as a constant value

* Adding models and constant values in this way could just as easily be done in a run script, but the advantage of defining a class is that it can be saved in a file (i.e. 'my_custom_phases') and reused in any project:
.. code-block:: Python

    from my_custom_phases import Oil
    oil = Oil(network=pn)

===============================================================================
Define Physics Objects for Each Geometry and Each Phase
===============================================================================

In the `previous tutorial <intermediate_usage>`_ we created two **Physics** object, one for each of the two **Geometry** objects used to handle the stratified layers.  In this tutorial, the internal pores and the boundary pores each have their own **Geometry**, but there are two **Phases**, which also each require a unique **Physics**:

.. code-block:: Python

    >>> phys_water_internal = op.Physics.GenericPhysics(network=pn, phase=water, geometry=geom)
    >>> phys_air_internal = op.Physics.GenericPhysics(network=pn, phase=air, geometry=geom)
    >>> phys_water_boundary = op.Physics.GenericPhysics(network=pn, phase=water, geometry=boun)
    >>> phys_air_boundary = op.Physics.GenericPhysics(network=pn, phase=air, geometry=boun)

* To reiterate, *one* **Physics** object is required for each **Geometry** *AND* each **Phase**, so the number can grow to become annoying very quickly  Some useful tips for easing this situation are given below.

-------------------------------------------------------------------------------
Create a Custom Pore-Scale Physics Model
-------------------------------------------------------------------------------

Perhaps the most distinguishing feature between pore-network modeling papers is the pore-scale physics models employed.  Accordingly, OpenPNM was designed to allow for easy customization in this regard, so that you can create your own models to augment or replace the ones included in the OpenPNM *models* libraries.  For demonstration, let's implement the capillary pressure model proposed by `Mason and Morrow in 1994 <http://dx.doi.org/10.1006/jcis.1994.1402>`_.  They studied the entry pressure of non-wetting fluid into a throat formed by spheres, and found that the converging-diverging geometry increased the capillary pressure required to penetrate the throat.  As a simple approximation they proposed :math:`P_c = -2 \sigma \cdot cos(2/3 \theta) / R_t`.

Pore-scale models are written as basic function definitions:

.. code-block:: Python
    :linenos:
    :caption: **Example of a Pore-Scale Model Definition**

    >>> def mason_model(network, phase, physics, f=0.6667, **kwargs):
    ...     Dt = network['throat.diameter']
    ...     theta=phase['throat.contact_angle']
    ...     sigma=phase['throat.surface_tension']
    ...     Pc = -4*sigma*sp.cos(f*sp.deg2rad(theta))/Dt
    ...     return Pc[network.throats(physics.name)]

Let's examine the components of above code:

* The function receives ``network``, ``phase`` objects as arguments.  Each of these provide access to the properties necessary for the calculation: ``'pore.diameter'`` values are retrieved via the ``network``, and the thermophysical properties are retrieved directly from the ``phase``.

* The ``f`` value is a scale factor that is applied to the contact angle.  Mason and Morrow suggested a value of 2/3 as a decent fit to the data, but we'll make this an adjustable parameter with 2/3 as the default.

* Note the ``pore.diameter`` is actually a **Geometry** property, but it is retrieved via the network using the data exchange rules outlined in the second tutorial, and explained fully in :ref:`data_storage`.

* All of the calculations are done for every throat in the network, but this pore-scale model is meant to be assigned to a single **Physics** object.  As such, the last line extracts values from the ``Pc`` array for the location of ``physics`` and returns just the subset.

* The actual values of the contact angle, surface tension, and throat diameter are NOT sent in as numerical arrays, but rather as dictionary keys to the arrays.  There is one very important reason for this: if arrays had been sent, then re-running the model would use the same arrays and hence not use any updated values.  By having access to dictionary keys, the model actually looks up the current values in each of the arrays whenever it is run.

* It would be a better practice to include the dictionary keys as arguments, such as ```contact_angle = 'throat.contact_angle'```.  This way the user could control where the contact angle could be stored on the **Phase** object.

Assuming this function is saved in a file called 'my_models.py' in the current working directory, this model can be used as:

.. code-block:: python

    from my_models import mason_model

-------------------------------------------------------------------------------
Copy Models Between Physics Objects
-------------------------------------------------------------------------------

As mentioned above, the need to specify a separate **Physics** object for each **Geometry** and **Phase** can become tedious.  It is possible to *copy* the pore-scale models assigned to one object onto another object.  First, let's assign the models we need to ``phys_water_internal``:

.. code-block:: python

    >>> phys_water_internal.models.add(propname='throat.capillary_pressure',
    ...                                model=mason_model)
    >>> mod = OpenPNM.Physics.hydraulic_conductance.hagan_poisseuille
    >>> phys_water_internal.models.add(propname='throat.hydraulic_conductance',
    ...                                model=mod)

Now make a copy of the ``models`` on ``phys_water_internal`` and apply it all the other **Physics** objects:

.. code-block:: python

    >>> mods = phys_water_internal.models.copy()
    >>> phys_water_boundary.models = mods
    >>> phys_air_internal.models = mods
    >>> phys_air_internal.models = mods

The only 'gotcha' with this approach is that each of the **Physics** objects must be *regenerated* in order to place numerical values for all the properties into the data arrays:

.. code-block:: python

    >>> phys_water_boundary.models.regenerate()
    >>> phys_air_internal.models.regenerate()
    >>> phys_air_internal.models.regenerate()

-------------------------------------------------------------------------------
Access Other Objects via the Network
-------------------------------------------------------------------------------

The above code used 3 lines to explicitly regenerate each **Physics** object, but an alternative and more efficient approach is possible.  When every object is created, it is 'registered' with the **Network** which is a required argument in the instantiation of every other object.  Any object can be looked-up by it's type using ``pn.geometries``, ``pn.phases``, or ``pn.physics``, which return a *dict* containing *key-value* pair of ``{object.name: object}``:

.. code-block:: python

    >>> pn.geometries
    {'internal': <OpenPNM.Geometry.__Stick_and_Ball__.Stick_and_Ball object at 0x7ef4e58>, 'boundary': <OpenPNM.Geometry.__GenericGeometry__.GenericGeometry object at 0x47734f8>}
    >>> pn.geometries.keys()  # Or obtain a list of object names using keys
    ['internal', 'boundary']

One handy use of this list is that is can be iterated over to perform an action on all objects in one line.  In this case running the ``regenerate`` method on all **Physics** objects can be accomplished with:

.. code-block:: python

    >>> temp = [item.regenerate for item in pn.physics.values()]

The ``values`` method of the *dict* class returns a list of the objects stored under each key.

-------------------------------------------------------------------------------
Adjust Pore-Scale Model Parameters
-------------------------------------------------------------------------------

The pore-scale models are stored in a **ModelsDict** object that is itself stored under the ``models`` attribute of each object.  This arrangement is somewhat convoluted, but it enables integrated storage of models on the object's wo which they apply.  The models on an object can be inspected with ``print(phys_water)``, which shows a list of all the pore-scale properties that are computed by a model, and some information about the model's *regeneration* mode.

Each model in the **ModelsDict** can be individually inspected by accessing it using the dictionary key corresponding to *pore-property* that it calculates, i.e. ``print(phys_water)['throat.capillary_pressure'])``.  This shows a list of all the parameters associated with that model.  It is possible to edit these parameters directly:

.. code-block:: python

    >>> phys_water.models['throat.capillary_pressure']['f']  # Inspect present value
    0.6666666666666666
    >>> phys_water.models['throat.capillary_pressure']['f'] = 0.75  # Change value

More details about the **ModelsDict** and **ModelWrapper** classes can be found in :ref:`models`.

===============================================================================
Perform Multiphase Transport Simulations
===============================================================================

-------------------------------------------------------------------------------
Use the Built-In Drainage Algorithm to Generate an Invading Phase Configuration
-------------------------------------------------------------------------------

.. code-block:: python

    >>> inv = op.Algorithms.Drainage(network=pn)
    >>> inv.setup(invading_phase=water, defending_phase=air)
    >>> inv.set_inlets(pores=pn.pores('top', 'bottom'))
    >>> inv.run()

* The inlet pores were set to both ``'top'`` and ``'bottom'`` using the ``pn.pores`` method.  The algorithm applies to the entire network so the mapping of network pores to the algorithm pores is 1-to-1.

* The ``run`` method automatically generates a list of 25 capillary pressure points to test, but you can also specify more pores, or which specific points to tests.  See the methods documentation for the details.

* Once the algorithm has been run, the resulting capillary pressure curve can be viewed with ``plot_drainage_curve``.  If you'd prefer a table of data for plotting in your software of choice you can use ``get_drainage_data`` which prints a table in the console.

-------------------------------------------------------------------------------
Set Pores and Throats to Invaded
-------------------------------------------------------------------------------

After running, the ``mip`` object possesses an array containing the pressure at which each pore and throat was invaded, stored as ``'pore.inv_Pc'`` and ``'throat.inv_Pc'``.  These arrays can be used to obtain a list of which pores and throats are invaded by water, using Boolean logic:

.. code-block:: python

    >>> Pi = inv['pore.inv_Pc'] < 10000
    >>> Ti = inv['throat.inv_Pc'] < 10000

The resulting Boolean masks can be used to manually adjust the hydraulic conductivity of pores and throats based on their phase occupancy.  The following lines set the water filled throats to near-zero air conductivity and vice-versa.

.. code-block:: python

    >>> phys_water['throat.hydraulic_conductance'][~Ti] = 1e-20
    >>> phys_air['throat.hydraulic_conductance'][Ti] = 1e-20

* The logic of these statements implicitly assumes that transport between two pores is only blocked if the throat is filled with the other phase, meaning that both pores could be filled and transport is still permitted.  Another option would be to set the transport to near-zero if *either* or *both* of the pores are filled as well.

* The above approach can get complicated if there are several **Geometry** objects, and it is also a bit laborious.  There is a pore-scale model for this under **Physics.models.multiphase** called ``conduit_conductance``.  The term conduit refers to the path between two pores that includes 1/2 of each pores plus the connecting throat.

-------------------------------------------------------------------------------
Calculate Relative Permeability of Each Phase
-------------------------------------------------------------------------------

We are now ready to calculate the relative permeability of the domain under partially flooded conditions.  Instantiate an **StokesFlow** object:

.. code-block:: python

    >>> water_flow = op.Algorithms.StokesFlow(network=pn, phase=water)
    >>> water_flow.set_boundary_conditions(pores=pn.pores('left'), bcvalue=200000, bctype='Dirichlet')
    >>> water_flow.set_boundary_conditions(pores=pn.pores('right'), bcvalue=100000, bctype='Dirichlet')
    >>> water_flow.run()
    >>> Q_partial = water_flow.rate(pores=pn.pores('right'))

The *relative* permeability is the ratio of the water flow through the partially water saturated media versus through fully water saturated media; hence we need to find the absolute permeability of water.  This can be accomplished by *regenerating* the ``phys_water`` object, which will recalculate the ``'throat.hydraulic_conductance'`` values and overwrite our manually entered near-zero values from the ``inv`` simulation using ``phys_water.models.regenerate()``.  We can then re-use the ``water_flow`` algorithm:

.. code-block:: python

    >>> water_flow.run()
    >>> Q_full = water_flow.rate(pores=pn.pores('right'))

And finally, the relative permeability can be found from:

.. code-block:: python

    >>> K_rel = Q_partial/Q_full

* The ratio of the flow rates gives the normalized relative permeability since all the domain size, viscosity and pressure differential terms cancel each other.

* To generate a full relative permeability curve the above logic would be placed inside a for loop, with each loop increasing the pressure threshold used to obtain the list of invaded throats (``Ti``).

* The saturation at each capillary pressure can be found be summing the pore and throat volume of all the invaded pores and throats using ``Vp = geom['pore.volume'][Pi]`` and ``Vt = geom['throat.volume'][Ti]``.

===============================================================================
Save the Simulation in a *PNM* File for Later Use
===============================================================================

OpenPNM includes a **Workspace** class that provides the type of functionality found on the *menu-bar* of a typical application GUI. Specifically, this enables *saving* and *loading* of all active networks, or individual objects.

To use these feature it is necessary to instantiate an instance:

.. code-block:: python

    >>> mgr = op.Base.Workspace()
    >>> mgr.save('filename.pnm')

Some of the more common functions of the **Workspace** are available via short-cuts under the main package, such that ``op.save`` is equivalent to calling ``mgr.save``.
