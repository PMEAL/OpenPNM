.. _data_storage:

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Data Storage
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Each OpenPNM Core object is a Python *dictionary* which is similar to a structured variable or *struct* in other languages.  This allows data to be stored and accessed by name, with a syntax like ``network['pore.diameter']``.  All pore and throat data are stored as vectors or lists, of either *Np* or *Nt* length representing the number of pores and throats on the object, respectively.  This means that each pore (or throat) has a number, and all properties for that pore (or throat) are stored in the array corresponding to that number.  Thus, the diameter for pore 15 is stored in the ``'pore.diameter'`` array in element 15, and the length of throat 32 is stored in the ``'throat.length'`` array at element 32.

This list-based approach is ideal when using the Numpy and Scipy libraries which are designed for elementwise, vectorized programming.  For instance, the volume of each throats can be found simultaneously using ``T_vol = 3.1415*(network['throat.radius']**2) * network['throat.length']``.  ``T_vol`` will be an *Nt*-long array of values, assuming ``'throat.length'`` and ``'throat.radius'`` were also *Nt*-long.

Several rules have been implemented to control the integrity of the data:

1. All array names must begin with either *'pore.'* or *'throat.'* which serves to identify the type of information they contain.
2. For the sake of consistency only arrays of length *Np* or *Nt* are allowed in the dictionary. Assigning a scalar value to a dictionary results in the creation of a full length vector, either *Np* or *Nt* long, depending on the name of the array..  This effectively applies the scalar value to all locations in the network.
3. Any Boolean data will be treated as a *label* while all other numerical data is treated as a *property*.  The difference between these is outlined below.

===============================================================================
Properties and Labels
===============================================================================
OpenPNM differentiates between two types of data for pores and throats: *properties* and *labels*.  The only difference between these is that *labels* are Boolean arrays (True / False), while *properties* are numerical data types.

-------------------------------------------------------------------------------
Properties, aka Data
-------------------------------------------------------------------------------
The physical details about pores and throats is referred to as *properties*, which includes information such as *pore volume* and *throat length*.  Properties are accessed using Python dictionary syntax to access the array of choice, then Numpy array indexing to access the pore or throat locations of choice:

.. code-block:: python

    >>> import OpenPNM
    >>> import scipy as sp
    >>> pn = OpenPNM.Network.Cubic(shape=[3, 3, 3])
    >>> pn['pore.coords'][1]
    array([ 0.5,  0.5,  1.5])

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
    ['pore.index', 'pore.coords', 'throat.conns']
    >>> pn.props('pore')
    ['pore.index', 'pore.coords']

You can also view a nicer list of ``props`` with ``print(pn.props())``.

-------------------------------------------------------------------------------
Labels
-------------------------------------------------------------------------------
Labels are a means of dynamically creating groups of pores and throats so they can be quickly accessed by the user.  For instance, is helpful to know which pores are on the *'top'* surface.  This label is automatically added by the *Cubic* network generator, so a list of all pores on the *'top'* can be retrieved by simply querying which pores possess the label *'top'* using the ``pores`` method:

.. code-block:: python

    >>> pn.pores('top')
    array([ 2,  5,  8, 11, 14, 17, 20, 23, 26])

The only distinction between *labels* and *properties* is that *labels* are Boolean masks of True/False.  Thus a ``True`` in element 10 of the array ``'pore.top'`` means that the label *'top'* has been applied to pore 10.  Adding and removing existing labels to pores and throats is simply a matter of setting the element to ``True`` or ``False``.  For instance, to remove the label *'top'* from pore 2:

.. code-block:: python

    >>> pn['pore.top'][2] = False
    >>> pn['pore.top']
    array([ 5,  8, 11, 14, 17, 20, 23, 26])
    >>> pn['pore.top'][2] = True  # Re-apply label to pore 2

Creating a new label array occurs automatically if a Boolean array is stored on an object:

.. code-block:: python

    >>> pn['pore.dummy_1'] = sp.rand(27) < 0.5

A complication arises if you have a list of pore numbers you wish to labels, such as [3, 4, 5].  You must first create the label array with all ``False`` values, *then* assign ``True`` to the desired locations:

.. code-block:: python

    >>> pn['pore.dummy_2'] = False  # Automatically assigns False to every pore
    >>> pn['pore.dummy_2'][[3, 4, 5]] = True
    >>> pn.pores('dummy_2')
    array([3, 4, 5])

The *label* functionality basically works by using Scipy's ``where`` method to return a list of locations where the array is ``True``:

.. code-block:: python

    >>> sp.where(pn['pore.dummy_2'])[0]
    array([3, 4, 5])

The ``pores`` and ``throats`` methods offer several useful enhancements to this approach.  For instance, several labels can be queried at once:

.. code-block:: python

    >>> pn.pores(['top', 'dummy_2'])
    array([ 2,  3,  4,  5,  8, 11, 14, 17, 20, 23, 26])

And there is also a ``mode`` argument which can be used to apply *set theory* logic to the returned list:

.. code-block:: python

    >>> pn.pores(['top', 'dummy_2'], mode='intersection')
    array([5])

This *set* logic basically retrieves a list of all pores with the label ``'top'`` and a second list of pores with the label ``dummy_2``, and returns the ``'intersection'`` of these lists, or only pores that appear in both lists.

The ``labels`` method can be used to obtain a list of all defined labels. This method optionally accepts a list of *pores* or *throats* as an argument and returns only the *labels* that have been applied to the specified locations.

.. code-block:: python

    >>> pn.labels()
    ['pore.all', 'pore.back', 'pore.bottom', 'pore.front', 'pore.internal', 'pore.left', 'pore.right', 'pore.top', 'throat.all']

This results can also be viewed with ``print(pn.labels())``.

.. note:: **The Importance of the 'all' Label**

   All objects are instantiated with a ``'pore.all'`` and ``'throat.all'`` label.  These arrays are essential to the framework since they are used to define how long the 'pore' and 'throat' data arrays must be.  In other words, the ``__setitem__`` method checks to make sure that any 'pore' array it receives has the same length as ``'pore.all'``.

-------------------------------------------------------------------------------
Counts and Indices
-------------------------------------------------------------------------------
One of the most common questions about a network is "how many pores and throats does it have?"  This can be answered very easily with the ``num_pores`` and ``num_throats`` methods.  Because these methods are used so often, there are also shortcuts: ``Np`` and ``Nt``.

.. code-block:: python

    >>> pn.num_pores()
    27

It is also possible to *count* only pores that have a certain label:

.. code-block:: python

    >>> pn.num_pores('top')
    9

These counting methods actually work by counting the number of ``True`` elements in the given *label* array.

Another highly used feature is to retrieve a list of pores or throats that have a certain label applied to them, which is of course is the entire purpose of the *labels* concept.  To receive a list of pores on the *'top'* of the **Network**:

.. code-block:: python

    >>> pn.pores('top')
    array([ 2,  5,  8, 11, 14, 17, 20, 23, 26])

The ``pores`` and ``throats`` methods both accept a *'mode'* argument that allows for *set-theory* logic to be applied to the query, such as returning 'unions' and 'intersections' of locations.

Often, one wants a list of *all** pore or throat indices on an object, so there are shortcut methods for this: ``Ps`` and ``Ts``.
