*******************************************************************************
Getting Started
*******************************************************************************

===============================================================================
Installation (?)
===============================================================================

Some installation instructions

-------------------------------------------------------------------------------
Preparing Your System
-------------------------------------------------------------------------------
OpenPNM utilizes Python 2.7 with Scipy 0.14 and Matplotlib.  The code is hosted on Github so it is also useful to have a Git client installed that can download the repository.  OpenPNM outputs many of its results in the vtk format which is conveniently visualized with Paraview.  

The installation of these requirements on specific platforms is described below.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Windows Instructions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
The simplest way to get Python, Scipy and Matplotlib is to download the `WinPython <http://code.google.com/p/winpython/>`_ package.  This package also comes with `Spyder <http://code.google.com/p/spyderlib/>`_, which provides an integrated development environment (IDE) that is very similar to Matlab, which an editor, command console, variable explorer and so on combined into the same window.  

Another useful tool is `Sourcetree <http://sourcetreeapp.com>`_ which provides a good interface to Github hosted repositories.  

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Apple Instructions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
TODO

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Linux Instructions
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
If you are using Linux, then you should ensure your Python is up to date, and that you have the Scipy and Matplot libraries installed.  

-------------------------------------------------------------------------------
Setting Up the Code
-------------------------------------------------------------------------------
The source code for OpenPNM is hosted on `Github <http://github.com/PMEAL/OpenPNM>`_.  You can download this code as a zip file or 'clone' the repository into your Git client (as as `Sourcetree <http://sourcetreeapp.com>`_).  These files can be located at a location of your preference, such as C:\user\documents\code\OpenPNM.  


===============================================================================
Example: Working with Scripts
===============================================================================
1.  Open Spyder
2.  In Spyder, select the 'Browse a Working Directory' button and locate your OpenPNM folder and select.  Also click the 'Set as console's working directory' button.
3.  Open a new blank .py file in the editor window if there is not one already.  

The first thing you must do is tell Spyder to import the OpenPNM code so you have access to the functions and methods, so add

.. code-block:: Python
   import OpenPNM








