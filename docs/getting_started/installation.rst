.. _installation:

================================================================================
Installation Instructions
================================================================================

--------------------------------------------------------------------------------
Requirements
--------------------------------------------------------------------------------

OpenPNM relies on the *Scipy Stack*, which includes scipy, numpy, matplotlib, pandas, and *many* other useful scientific packages.  These packages are extremely tedious to install from source (i.e. doing ``pip install numpy`` is a terrible idea).  It is highly recommended to first download and install the `anaconda distribution <https://www.anaconda.com/download>`_ for your system.  Be sure to get the latest Python 3 version.  Once this is installed, you will then have a fully functioning and highly optimized compiled versions of numpy, scipy and the other great packages, and you'll be ready to use OpenPNM.

--------------------------------------------------------------------------------
Installing from PyPI
--------------------------------------------------------------------------------

OpenPNM can be installed using pip, as follows, which downloads the latest stable release from the Python Package Index:

.. code-block:: python

    pip install openpnm

This will install OpenPNM into your Python environment.

To upgrade your OpenPNM to a newer version, use ``pip install --upgrade openpnm``.

Alternatively, you might wish to use the latest development version, and perhaps to edit the source code yourself.  In this case, you should download or clone the repository from Github, and save it to a preferred location on your local machine.  Once this is done you can install from the git repository using:

.. code-block:: python

    pip install -e "path/to/downloaded/git/repo"

The -e commands make the installation editable, so that any changes made the source files will be available the next time OpenPNM is imported.  This is useful for switch branches in the Git repo, or for adding new algorithms to a custom branch.

--------------------------------------------------------------------------------
Other Requirements
--------------------------------------------------------------------------------
It is also suggested to download `Paraview <http://www.paraview.org/>`_ for visualizing the networks produced by OpenPNM.
