.. _installation:

===============================================================================
Installation Instructions
===============================================================================

Before proceeding with the installation, be sure to read the *Requirements* section below.

Installing OpenPNM is done by running the following on your command line:

.. code-block:: python

    pip install openpnm

This will install OpenPNM into your Python environment.  To use OpenPNM, open a Python console and type:

>>> import OpenPNM

To upgrade your OpenPNM to a newer version, use ``pip install --upgrade openpnm``.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Requirements
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
OpenPNM relies on the *Scipy Stack*, which include scipy, numpy, matplotlib, pandas, and *many* other useful scientific packages.  These packages are extremely tedious to install from source (i.e. doing ``pip install numpy`` is a terrible idea).  It is highly recommended to first download and install the [anaconda distribution](https://www.anaconda.com/download) for your system.  Be sure to get the Python 3+ version.  Once this is installed, you will then have fully functioning and highly optimized compiled versions of numpy, scipy and the other great packages.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Other Requirements
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
It is also suggested to download `Paraview <http://www.paraview.org/>`_ for visualizing the networks produced by OpenPNM.
