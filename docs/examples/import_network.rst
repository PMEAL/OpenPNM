.. _network_import:

===============================================================================
Importing Networks Created using External Code
===============================================================================

OpenPNM supports importing network data in several formats:
	* The recommended format is 'csv', which is human readable, editable in spreadsheet programs, and supported by almost all software since it's so simple.
	* The 'mat' format is supported for importing from Matlab, which is a common environment for many researchers.
	* OpenPNM uses Paraview to visualize output via the 'vtp' format, so we also support importing from these files. 'vtp' is a well defined and common format so it's also possible to export to 'vtp' from numerous programs.
	* Finally, OpenPNM aims to support specialized formats used by specific software and research groups.  Including:
		- The 'Statoil' format that is used by Martin Blunt's Maximal Ball network extraction code.
		- The 'yaml' format is specifically intended for importing from the NetworkX package.

.. note:: Import File Format Specifications

    The actual internal formatting of each file type is specified in the documentation for each of the IO classes.  These are located under OpenPNM.Utilities.IO.

-------------------------------------------------------------------------------
Using the Import Network Class
-------------------------------------------------------------------------------

There are two ways to import data onto an OpenPNM **Network**.  For the first approach, start by instantiating an *Import* class:

.. code-block:: python

    import OpenPNM
    import scipy as sp
    pn = OpenPNM.Network.Import()

You'll notice that ``pn`` has all the usual methods and functions of a **Network** object, but also has some ``from_\*`` methods.  These are used to import data 'from' the various supported file formats.  These can be used to incrementally add data from different sources onto the same network, so if you import a Statoil network, but have also performed some additional calculations in Matlab, you can use ``from_mat`` to add this extra data to ``pn`` without overwriting the existing data (assuming the ``mode`` argument is set to 'add' not 'overwrite').

For the first part of this tutorial, we'll import a netowrk saved in CSV format, which can be `found here <https://db.tt/yOSJQcds>`_.  Download this file to your working directory, and load it into ``pn`` using:

.. code-block:: python

    pn.from_csv(filename='test_load_csv')  # doctest: +SKIP

You can inspect this network using the usual methods, such as ``num_pores`` and ``print(pn)``.  This particular CSV file was filled with a full set of geometrical properties, so the network is nearly complete.  You can of course add and remove data and labels as needed, and also add pore-scale models as well.

.. code-block:: python

    del pn['throat.surface_area']  # doctest: +SKIP
    pn['throat.surface_area'] = sp.constants.pi/4*pn['throat.length']*pn['throat.diameter']**2

-------------------------------------------------------------------------------
Using the IO Classes Directly
-------------------------------------------------------------------------------

The *Import* network class described above is a 'wrapper' class for the actual import/export classes which are stored in OpenPNM.Utilities.IO.  The second method for importing network data uses these classes directly.

.. code-block:: python

    pn = OpenPNM.Utilities.IO.CSV.load(filename='test_load_csv')

		
