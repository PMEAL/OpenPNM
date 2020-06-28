.. _installation:

Installing OpenPNM
##################

OpenPNM requires the *Scipy Stack* (Numpy, Scipy, Matplotlib, etc), which is most conveniently obtained by installing the `Anaconda Distribution <https://conda.io/docs/user-guide/install/download.html>`_, so first download and install this if you don't have it already.

OpenPNM can be installed from the `Python Package Index https://pypi.org/project/openpnm/`_ using:

.. code-block::

   pip install openpnm

Or the source code can be downloaded from `Github <https://github.com/pmeal/OpenPNM/>`_ and installed by running:

.. code-block::
   pip install -e 'path/to/downloaded/files'

The advantage to installing from the source code is that you can edit the files and have access to your changes each time you import *OpenPNM*.

Finally, OpenPNM is now available as a `conda package <https://anaconda.org/conda-forge/openpnm>`_, so you can also do:

.. code-block::

   conda install -c conda-forge openpnm
