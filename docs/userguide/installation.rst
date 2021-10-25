.. _installation:


Installation
------------

OpenPNM requires the *Scipy Stack* (Numpy, Scipy, Matplotlib, etc), which is most conveniently obtained by installing the `Anaconda Distribution <https://conda.io/docs/user-guide/install/download.html>`_, so first download and install this if you don't have it already.


Preferred method
----------------
The preferred way of installing OpenPNM is via Anaconda Cloud using:

.. code-block::

    conda install -c conda-forge openpnm

Alternative method
------------------
OpenPNM can also be installed from the `Python Package Index <https://pypi.org/project/openpnm/>`_ using:

.. code-block::

    pip install openpnm

However, we don't recommend installing using pip since ``pypardiso``, which is a blazing fast direct solver, is not available for Windows users who use Python 3.7+.

For developers
--------------
For developers who intend to change the source code or contribute to OpenPNM, the source code can be downloaded from Github and installed by running:

.. code-block::

    pip install -e 'path/to/downloaded/files'

The advantage to installing from the source code is that you can edit the files and have access to your changes each time you import OpenPNM.
