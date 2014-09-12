.. _installation:

===============================================================================
Installation Instructions 
===============================================================================
Installing OpenPNM is done by running the following on your command line:

.. code-block:: python

    pip install openpnm

This will install OpenPNM into your Python environment.  To use OpenPNM, open a Python console and type:

>>> import OpenPNM

You can check your version of OpenPNM with:

>>> OpenPNM.__version__
'1.0.0'

To upgrade your OpenPNM to a newer version, use ``pip install --upgrade openpnm``.

It is also possible to download the source code directly from Github and work from that.  This is not recommended unless you are planning to do 'development' work on the framework.  The pip install approach places the source code in the Python directory and out of harms way.  

===============================================================================
Installing a Python Environment
===============================================================================
The above installation instructions assume that you have Python installed on your system, and the ``pip`` installer.  If you do not, then refer to the sections below for your system.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Windows
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
To use OpenPNM it is necessary to have Python installed on your system.  The simplest way to get Python and all the helpful add-ons. is to download the `WinPython <http://code.google.com/p/winpython/>`_ package.  This package also comes with `Spyder <http://code.google.com/p/spyderlib/>`_, which provides an integrated development environment (IDE) that is very similar to Matlab, with an editor, command console, variable explorer and so on combined into the same window.  

Once WinPython is install, navigate to the directory where you chose to install it and open spyder.exe.  From the 'tools' method, select 'open command prompt'.  This will open a special version of the Window command prompt that is aware of the Python installation.  Here you will use ``pip install openpnm`` to have OpenPNM installed on your machine.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Apple Instructions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The *Conda* package by `Continuum Analytics <http://continuum.io/downloads#all?`_ is probably the best option for Mac users.  This will download Python and all the usual scientific packages, as well as Spyder which the OpenPNM developers highly recommend.  

Once Conda is installed, you can navigate to the Conda directory and open the Conda command prompt application.  Conda is a little bit more convoluted that WinPython because they offer a special version of ``pip`` called ``conda``.  OpenPNM is not registed with the ``conda`` repository so you must use ``pip install openpnm``.  This will work, but we've found that it does not play nice with virtual environments.  As long as you are using Conda without being too fancy things should work fine.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Linux Instructions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
If you are using Linux, then you probably don't need instructions :-)

===============================================================================
Other Requirements
===============================================================================
It is also suggested to download `Paraview <http://www.paraview.org/>`_ for visualizing the networks produced by OpenPNM.