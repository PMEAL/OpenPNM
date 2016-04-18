.. _advanced_tutorial3

###############################################################################
Tutorial 3 of 3: Advanced Topics and Usage
###############################################################################

In this final tutorial on the use of OpenPNM, the focus will be on creating custom pore-scale models and custom classes, and manipulating topology

.. contents:: Topics Covered in this Tutorial

**Learning Outcomes**

#. Use different methods to add boundary pores to a network
#. Manipulate network topology by adding and removing pores and throats
#. Explore the ModelsDict design, including copying models between objects, and changing model parameters
#. Write a custom pore-scale model
#. Retrieve objects of a certain type associated with the current network
#. Combine multiple algorithms to predict relative permeability
#. Use the workspace manager to save and load, clone and purge simulations

===============================================================================
Build and Manipulate Network Topology
===============================================================================

For the present tutorial, we'll keep the topology simple to help keep the focus on other aspects of OpenPNM.

.. code-block:: python

    >>> import scipy as sp
    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[10, 10, 10], spacing=0.00006)

-------------------------------------------------------------------------------
Adding Boundary Pores
-------------------------------------------------------------------------------

When performing transport simulations it is often useful to have 'boundary' pores attached to the surface(s) of the network where boundary conditions can be applied.  The **Cubic** class has two methods for doing this: ``add_boundaries`` and ``add_boundary_pores``.  The first method adds boundaries to ALL six faces of the network and offsets them from the network by 1/2 of the value provided as the network ``spacing``.  The second method provides total control over which boundaries are created and where they are positioned, but it more cumbersome to use.  Let's explore these two options:

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

Let's remove all these newly created boundary pores.  When they are created these pores are all automatically labeled with a label such as ``'top_boundary'``, so we can select all boundary pores using:

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

    >>> Ts = sp.rand(pn.Nt) < 0.2  # Create a mask with ~20% of throats labeled True
    >>> pn.trim(throats=Ts)  # Use mask to indicate which throats to trim

When the ``trim`` function is called, it automatically checks the health of the network afterwards, so logger messages might appear on the command line if problems were found.  **Networks** have a ``check_network_health`` method that performs the same checks, and returns a **HealthDict** containing the results of the checks:

.. code-block:: python

    >>> a = pn.check_network_health()
    >>> pn.trim(pores=a['trim_pores'])

The **HealthDict** contains several lists including things like duplicate throats and isolated pores, but also a suggestion of which pores to trim to return the network to a healthy state.

===============================================================================
Define Geometry Objects
===============================================================================

Since we've added boundary pores to the network we need to the treat them a little bit differently.  Specifically, they should have no volume or length (as they are not physically representative of real pores).  To do this, we create two separate **Geometry** objects, one for internal pores and one for the boundaries:

.. code-block:: python

    >>> Ps = pn.pores('*boundary', mode='not')
    >>> geom = op.Geometry.Stick_and_Ball(network=pn, pores=Ps, throats=pn.Ts)
    >>> Ps = pn.pores('*boundary')
    >>> boun = op.Geometry.GenericGeometry(network=pn, pores=Ps)

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

    >>> water = OpenPNM.Phases.Water(network=pn)
    >>> air = OpenPNM.Phases.Air(network=pn)

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

.. code-block:: python

    >>> phys_water_internal = OpenPNM.Physics.GenericPhysics(network=pn, phase=water, geometry=geom)
    >>> phys_air_internal = OpenPNM.Physics.GenericPhysics(network=pn, phase=air, geometry=geom)
    >>> phys_water_boundary = OpenPNM.Physics.GenericPhysics(network=pn, phase=water, geometry=boun)
    >>> phys_air_boundary = OpenPNM.Physics.GenericPhysics(network=pn, phase=air, geometry=boun)

* To reiterate, *one* **Physics** object is required for each **Geometry** *AND* each **Phase**, so the number can grow to become annoying very quickly  Some useful tips for easing this situation are given below.

-------------------------------------------------------------------------------
Create a Custom Pore-Scale Physics Model
-------------------------------------------------------------------------------

Perhaps the most distinguishing feature between pore-network modeling papers is the pore-scale physics models employed.  OpenPNM was designed to allow for easy customization in this regard, so that you can create your own models to augment or replace the ones included in the OpenPNM *models* libraries.  For demonstration, let's implement the capillary pressure model proposed by `Mason and Morrow in 1994 <http://dx.doi.org/10.1006/jcis.1994.1402>`_.  They studied the entry pressure of non-wetting fluid into a throat formed by spheres, and found that the converging-diverging geometry increased the capillary pressure required to penetrate the throat.  As a simple approximation they proposed :math:`P_c = -2 \sigma \cdot cos(2/3 \theta) / R_t`.

Pore-scale models are written as basic function definitions:

.. code-block:: python
    :linenos:
    :caption: **Example of a Pore-Scale Model Definition**

    def mason_model(network, phase, physics,
                    contact_angle='throat.contact_angle',
                    surface_tension='throat.surface_tension',
                    diameter='throat.diameter',
                    **kwargs):
        Dt = network[diameter]
        theta=phase[contact_angle]
        sigma=phase[surface_tension]
        Pc = -4*sigma*sp.cos(2/3*sp.deg2rad(theta))/Dt
        return Pc[network.throats(physics.name)]

Let's examine the components of above code:

* The function receives ``network``, ``phase`` objects as arguments.  Each of these provide access to the properties necessary for the calculation.  The ``'pore.diameter'`` values are retrieved via the ``network``, and the thermophysical properties are retrieved directly from the ``phase``.

* Note the ``pore.diameter`` is actually a **Geometry** property, but it is retrieved via the network using the data exchange rules outlined in the second tutorial, and explained fully in :ref:`data_storage`.

* All of the calculations are done for every throat in the network, but this pore-scale model is meant to be assigned to a single **Physics** object.  As such, the last line extracts value from the ``Pc`` array for the location of ``physics`` and returns just the subset.

* The actual values of the contact angle, surface tension, and throat diameter are NOT sent in as numerical arrays, but rather as dictionary keys to the arrays.  There is one very important reason for this: if arrays had been sent, then re-running the model would use the same arrays and hence not use any updated values.  By having access to dictionary keys, the model actually looks up the current values in each of the arrays whenever it is run.

Assuming this function is saved in a file called 'my_models.py' in the current working directory, this model can be used as:

.. code-block:: python

    from my_models import mason_model

-------------------------------------------------------------------------------
Copy Models between Physics Objects
-------------------------------------------------------------------------------

As mentioned above, the need to specify a separate **Physics** object for each **Geometry** and **Phase** can become tedious.  It is possible to *copy* the pore-scale models assigned to one object onto another object.  First, let's assign the models we need to ``phys_water_internal``:

.. code-block:: python

    >>> phys_water_internal(propname='throat.capillary_pressure',
    ...                     model=mason_model)
    >>> phys_water_internal(propname='throat.hydraulic_conductance',
    ...                     model=OpenPNM.Physics.hydraulic_conductance.hagan_poisseuille)

Now we can make a copy of the ``models`` on ``phys_water_internal`` and apply it all the other **Physics** objects:

.. code-block:: python

    >>> mods = phys_water_internal.models.copy()
    >>> phys_water_boundary.models = mods
    >>> phys_air_internal.models = mods
    >>> phys_air_internal.models = mods

-------------------------------------------------------------------------------
Access Objects of a Certain Type via the Network
-------------------------------------------------------------------------------

The only 'gotcha' with the above approach is that each of the **Physics** objects must be *regenerated* in order to place numerical values for all the properties into the data arrays.  It would be possible to add several new lines to your run script that calls the ``regenerate`` method for each new object.  A more efficient approach is to utilize the fact that when every object is created, it is 'registered' with the **Network**.  Any object can be looked-up by it's type using ``pn.geometries()``, ``pn.phases()``, or ``pn.physics()``.

.. code-block:: python

    >>> temp = [item.regenerate for item in pn.physics()]


-------------------------------------------------------------------------------
Adjust Pore-Scale Model Parameters
-------------------------------------------------------------------------------




===============================================================================
Perform Multiphase Transport Simulations
===============================================================================



-------------------------------------------------------------------------------
Use the Built-In Drainage Algorithm to Generate an Invading Phase Configuration
-------------------------------------------------------------------------------



-------------------------------------------------------------------------------
Set Pores and Throats to Invaded
-------------------------------------------------------------------------------



-------------------------------------------------------------------------------
Calculate Relative Permeability of Each Phase
-------------------------------------------------------------------------------

blah
