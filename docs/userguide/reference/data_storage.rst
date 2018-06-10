.. _data_storage:

###############################################################################
Data Storage
###############################################################################

Each OpenPNM Core object is a Python *dictionary* which is similar to a structured variable or *struct* in other languages.  This allows data to be stored and accessed by name, with a syntax like ``network['pore.diameter']``.   Inside each *dict* are stored numerous arrays containing pore or throat data corresponding to the *key* (i.e. ``'pore.diameter'`` values).

All pore and throat data are stored in arrays of either *Np* or *Nt* length representing the number of pores and throats on the object, respectively.  This means that each pore (or throat) has a number that is implicitly indicated by it's location in the arrays.  All properties for pore *i* or throat *j* are stored in the array at the element *i* or *j*.  Thus, the diameter for pore 15 is stored in the ``'pore.diameter'`` array in element 15, and the length of throat 32 is stored in the ``'throat.length'`` array at element 32.  This array-based approach is ideal when using the Numpy and Scipy libraries which are designed for elementwise, vectorized programming.  For instance, the volume of each throats can be found simultaneously using ``T_vol = 3.1415*(network['throat.radius']**2) * network['throat.length']``.  ``T_vol`` will be an *Nt*-long array of values, assuming ``'throat.length'`` and ``'throat.radius'`` were also *Nt*-long.

Several rules have been implemented to control the integrity of the data:

#. All array names must begin with either *'pore.'* or *'throat.'* which serves to identify the type of information they contain.
#. For the sake of consistency only arrays of length *Np* or *Nt* are allowed in the dictionary. Assigning a scalar value to a dictionary results in the creation of a full length vector, either *Np* or *Nt* long, depending on the name of the array..  This effectively applies the scalar value to all locations in the network.
#. Any Boolean data will be treated as a *label* while all other numerical data is treated as a *property*.  The difference between these is outlined below.

===============================================================================
Properties (aka Data)
===============================================================================

The physical details about pores and throats are referred to as *properties*, which includes information such as *pore volume* and *throat length*.  Properties are accessed using Python dictionary syntax to access the array of choice, then Numpy array indexing to access the pore or throat locations of choice:

.. code-block:: python

    >>> import OpenPNM
    >>> import scipy as sp
    >>> pn = OpenPNM.Network.Cubic(shape=[3, 3, 3])
    >>> pn['pore.coords'][1]
    array([0.5, 0.5, 1.5])

Note that ``pn['pore.coords']`` retrieves the Numpy array from the dictionary, while the ``[1]`` retrieves the value in element 1 of the Numpy array.

Writing data is straightforward:

.. code-block:: python

    >>> pn['pore.foo'] = 1.0
    >>> pn['pore.foo'][5]
    1.0
    >>> pn['pore.foo'][6] = 2.0
    >>> pn['pore.foo'][6]
    2.0
    >>> pn['pore.foo'][5]
    1.0

The above lines illustrate how a scalar value is converted to a vector (*Np*-long), and how specific pore values can be assigned.  It is also possible to assign an entire array in one step:

.. code-block:: python

    >>> pn['pore.bar'] = sp.rand(27)  # pn has 27 pores (3*3*3)

Attempts to write an array of the wrong size will result in an error:

.. code-block:: python

    >>> pn['pore.baz'] = [2, 3, 4]

To quickly see a complete list *properties* on an object use the ``props`` method.  You can specify whether only *pore* or *throat* properties should be returned, but the default is both:

.. code-block:: python

    >>> pn.props()
    ['pore.bar', 'pore.coords', 'pore.foo', 'pore.index', 'throat.conns']
    >>> pn.props('throat')
    ['throat.conns']

You can also view a nicely formatted list of ``props`` with ``print(pn.props())``.

===============================================================================
Labels
===============================================================================
Labels are a means of dynamically creating groups of pores and throats so they can be quickly accessed by the user.  For instance, is helpful to know which pores are on the *'top'* surface.  This label is automatically added by the *Cubic* network generator, so a list of all pores on the *'top'* can be retrieved by simply querying which pores possess the label *'top'* using the ``pores`` method:

.. code-block:: python

    >>> pn.pores('top')
    array([ 2,  5,  8, 11, 14, 17, 20, 23, 26])

The only distinction between *labels* and *properties* is that *labels* are Boolean masks of True/False.  Thus a ``True`` in element 10 of the array ``'pore.top'`` means that the label *'top'* has been applied to pore 10.  Adding and removing existing labels to pores and throats is simply a matter of setting the element to ``True`` or ``False``.  For instance, to remove the label *'top'* from pore 2:

.. code-block:: python

    >>> pn['pore.top'][2] = False
    >>> list(sp.where(pn['pore.top'])[0])
    [5, 8, 11, 14, 17, 20, 23, 26]
    >>> pn['pore.top'][2] = True  # Re-apply label to pore 2

Creating a new label array occurs automatically if a Boolean array is stored on an object:

.. code-block:: python

    >>> pn['pore.dummy_1'] = sp.rand(27) < 0.5

A complication arises if you have a list of pore numbers you wish to label, such as [3, 4, 5].  You must first create the label array with all ``False`` values, *then* assign ``True`` to the desired locations:

.. code-block:: python

    >>> pn['pore.dummy_2'] = False  # Automatically assigns False to every pore
    >>> pn['pore.dummy_2'][[3, 4, 5]] = True
    >>> list(pn.pores('dummy_2'))
    [3, 4, 5]

The *label* functionality uses Scipy's ``where`` method to return a list of locations where the array is ``True``:

.. code-block:: python

    >>> list(sp.where(pn['pore.dummy_2'])[0])
    [3, 4, 5]

The ``pores`` and ``throats`` methods offer several useful enhancements to this approach.  For instance, several labels can be queried at once:

.. code-block:: python

    >>> list(pn.pores(['top', 'dummy_2']))
    [2, 3, 4, 5, 8, 11, 14, 17, 20, 23, 26]

And there is also a ``mode`` argument which can be used to apply *set theory* logic to the returned list:

.. code-block:: python

    >>> list(pn.pores(['top', 'dummy_2'], mode='intersection'))
    [5]

This *set* logic basically retrieves a list of all pores with the label ``'top'`` and a second list of pores with the label ``dummy_2``, and returns the ``'intersection'`` of these lists, or only pores that appear in both lists.

The ``labels`` method can be used to obtain a list of all defined labels. This method optionally accepts a list of *pores* or *throats* as an argument and returns only the *labels* that have been applied to the specified locations.

.. code-block:: python

    >>> pn.labels()
    ['pore.all', 'pore.back', 'pore.bottom', 'pore.dummy_1', 'pore.dummy_2', 'pore.front', 'pore.internal', 'pore.left', 'pore.right', 'pore.top', 'throat.all']

This results can also be viewed with ``print(pn.labels())``.

.. note:: **The Importance of the 'all' Label**

   All objects are instantiated with a ``'pore.all'`` and ``'throat.all'`` label.  These arrays are essential to the framework since they are used to define how long the 'pore' and 'throat' data arrays must be.  In other words, the ``__setitem__`` method checks to make sure that any 'pore' array it receives has the same length as ``'pore.all'``.

===============================================================================
Counts and Indices
===============================================================================
One of the most common questions about a network is "*how many pores and throats does it have?*"  This can be answered  easily with the ``num_pores`` and ``num_throats`` methods.  Because these methods are used so often, there are also shortcuts: ``Np`` and ``Nt``.

.. code-block:: python

    >>> pn.num_pores()
    27
    >>> pn.Np
    27

It is also possible to *count* only pores that have a certain label:

.. code-block:: python

    >>> pn.num_pores('top')
    9

These counting methods actually work by counting the number of ``True`` elements in the given *label* array.

Another highly used feature is to retrieve a list of pores or throats that have a certain label applied to them, which is of course is the entire purpose of the *labels* concept.  To receive a list of pores on the *'top'* of the **Network**:

.. code-block:: python

    >>> list(pn.pores('top'))
    [2, 5, 8, 11, 14, 17, 20, 23, 26]

The ``pores`` and ``throats`` methods both accept a *'mode'* argument that allows for *set-theory* logic to be applied to the query, such as returning 'unions' and 'intersections' of locations.

Often, one wants a list of *all** pore or throat indices on an object, so there are shortcut methods for this: ``Ps`` and ``Ts``.

It is also possible to filter a list of pores or throats according to their labels using ``filter_by_label``:

.. code-block:: python

    >>> Ps = pn.pores('top')
    >>> list(Ps)
    [2, 5, 8, 11, 14, 17, 20, 23, 26]
    >>> list(pn.filter_by_label(pores=Ps, labels='left'))
    [2, 11, 20]

The ``filter_by_label`` method also accepts a ``mode`` argument that applies additional filtering to the returned list using *set-theory*-type logic.  In this case, the method will find sets of pores or throats that satisfies each given label, then determines the *union*, *intersection*, or *difference* of the given sets.

===============================================================================
Data Exchange Between Objects
===============================================================================

One of the features in OpenPNM is the ability to model heterogeneous materials by apply different pore-scale models to different regions.  This is done by (a) creating a unique **Geometry** object for each region (i.e. small pores vs big pores) and (b) creating unique **Physics** object for each region as well (i.e. Knudsen diffusion vs Fickian diffusion).  One consequence of this segregation of properties is that a *single* array containing values for all locations in the domain cannot be directly obtained.  It is possible to manually piece together values from different regions, but this is cumbersome.  OpenPNM offers a shortcut for this, by making it possible to query **Geometry** properties via the **Network** object, and **Physics** properties from the associated **Phase** object:

::

    Documentation not finished yet
