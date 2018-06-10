.. _data_IO:

###############################################################################
Importing and Exporting Data
###############################################################################
OpenPNM has functions for importing and exporting data in various formats that are suitable for exchanging data

.. note::

    Exporting data to these file formats is **NOT** the recommended way to save and load OpenPNM Simulations, as lots of key information is lost (only numerical values are exported or imported).  To save and load simulations, use the methods available on the **Workspace** as described in :ref:`workspace`.

===============================================================================
Exporting Data
===============================================================================
OpenPNM allows for exporting the data to several formats:

#. CSV (Comma separated values) is the recommended format as it is widely used by almost all other software tools
#. MAT (Matlab file) is supported for the obvious reason that Matlab is a very popular choice for post-processing
#. VTK (Visualization Toolkit) is main format for exporting results to a visualization software (i.e. Paraview)

There are several ways to import and export data.  All the import and export classes are stored under ``OpenPNM.Utilities.IO``, but there is also ``import_data`` and ``export_data`` methods available in the top level of the project's namespace for convenience (i.e. ``OpenPNM.export_data``). The **Workspace** object also possess ``import_data`` and ``export_data`` methods.  All these approaches utilize the classes stored in the **Utilities.IO** module.

-------------------------------------------------------------------------------
Comma Separated Variables (or is it Values?)
-------------------------------------------------------------------------------
CSV files are the recommended format in OpenPNM due to their simplicity and wide interoperability with virtually all other software.  The list-type data storage scheme used in OpenPNM also fits very well in the CSV column-based format.

Exporting data is accomplished by:

.. code-block:: python


Exporting data is accomplished by:

.. code-block:: python

    >>> import OpenPNM as op
    >>> pn = op.Network.Cubic(shape=[5, 5, 5])
    >>> geom = op.Geometry.Stick_and_Ball(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> op.Utilities.IO.CSV.save(network=pn, filename='test_file.csv')

The resulting *csv* file contains all the data on the network ``pn``, but also from ``geom``.  The ability of **Networks** to retrieve data from its associated **Geometries**


-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Mat-File
-------------------------------------------------------------------------------

-------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Visualization Toolkit (VTK)
-------------------------------------------------------------------------------


===============================================================================
Importing Data
===============================================================================
OpenPNM supports more import formats than the export formats listed above. This is to accommodate the diverse range of external packages used to produce networks.  New import formats are added as the need arises.

# CSV, MAT and VTK are able to import data exported by OpenPNM as well as any other similarly formatted data
# Statoil is available specifically to import networks extracted using the maximal ball code of the Blunt group at ICL
# NetworkX is able to import networks from the NetworkX python package
