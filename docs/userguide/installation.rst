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

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Installing a Environment of the Python SciPy Stack 
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The above installation instructions assume that you have the Python SciPy Stack (Python including including packages like NumPy, SciPy, matplotlib) installed on your system, and the ``pip`` installer.  If you do not, then refer to the sections below for your system.  OpenPNM is designed to run with Python 3. 

-------------------------------------------------------------------------------
Windows Instructions
-------------------------------------------------------------------------------
The simplest way to get Python and all the necessary Python packages for scientific computing in Windows is to download the `WinPython <http://code.google.com/p/winpython/>`_ package.  This package also comes with `Spyder <http://code.google.com/p/spyderlib/>`_, which provides an integrated development environment (IDE) that is very similar to Matlab, with an editor, command console, variable explorer and so on combined into the same window.  

Once WinPython is installed, navigate to the directory where you chose to install it and open spyder.exe.  From the 'tools' method, select 'open command prompt'.  This will open a special version of the Window command prompt that is aware of the Python installation.  Here you will use ``pip install openpnm`` to have OpenPNM installed on your machine.  

As an alternative, you can also use the *Anaconda* distribution by `Continuum Analytics <http://continuum.io/downloads#all?>`_ for Windows (cf. below for more details).

-------------------------------------------------------------------------------
Apple Instructions
-------------------------------------------------------------------------------
The *Anaconda* distribution by `Continuum Analytics <http://continuum.io/downloads#all?>`_ is probably the best option for Mac users.  Download the installer and follow the install instructions as stated. This will provide you with Python and all the usual scientific packages, as well as Spyder which the OpenPNM developers highly recommend.  

Once Anaconda is installed, you can start the Anaconda launcher by double clicking ``Launcher.app`` in your ``~/anaconda`` directory, by which you can start a Spyder session, or also launch Spyder directly from the command line by running ``spyder``.  Anaconda comes with its own package management ``conda``, which you can run from the command line to install specific packages or to setup virtual environments. OpenPNM is not registed with the ``conda`` repository so you must use ``pip install openpnm``.  Using ``conda`` you can also update your packages, either updating all by running ``conda update anaconda``, which updates all of your packages, or update individual packages, e.g. updating the Python package itself by running ``conda upate python``.
Anaconda offers the option of running virtual environments next to each other. This allows to easily switch between different versions of Python packages, i.e. having one Python 2.7 enviroment and one Python 3.3 environment. This is possible by running ``conda create -n py3k python=3.4 numpy=1.8.1 scipy=0.14 anaconda`` from a command line, where you can explicitly specify which package versions should be installed, and now created a new environment called ``py3k``. In order to activate this virtual environment, run ``source activate py3k`` from a command line, and the run ``spyder``, which will give you a Spyder session using these specific packages. To run OpenPNM in this environment, you need to first activate the environment with ``source activate py3k`` and then run ``pip install openpnm``, so that OpenPNM is installed in this specific environment.

-------------------------------------------------------------------------------
Linux Instructions
-------------------------------------------------------------------------------
If you are using Linux, we also recommend the *Anaconda* distribution by `Continuum Analytics <http://continuum.io/downloads#all?>`. Download the installer and follow the install instructions as stated. This will provide you with Python and all the usual scientific packages, as well as Spyder which the OpenPNM developers highly recommend.  The process of using it is basically identical to Mac environment, so you can look at the instructions for Apple users stated above. The only difference is that the ``Launcher.app`` is adapted to Linux specifics, but also called ``Launcher``. Your Anaconda distribution is also installed in ``~/anaconda`` and you can use the ``conda`` commands as stated above to set up environments or update packages. 

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Other Requirements
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
It is also suggested to download `Paraview <http://www.paraview.org/>`_ for visualizing the networks produced by OpenPNM.