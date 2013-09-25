*******************************************************************************
Getting Started
*******************************************************************************

===============================================================================
Example: Working with Scripts
===============================================================================
1.  Open Spyder
2.  In Spyder, select the 'Browse a Working Directory' button and locate your OpenPNM folder and select.  Also click the 'Set as console's working directory' button.
3.  Open a new blank .py file in the editor window if there is not one already.  

The first thing you must do is tell Spyder to import the OpenPNM code so you have access to the functions and methods, so add

.. code-block:: Python
   import OpenPNM

This imports the OpenPNM package with all it's algorithms and data storage functions.

Next, you'll want to generate a network.  This is accomplished by first creating a generator, which we'll call gn1:

.. code-block:: Python
   gn1 = OpenPNM.Generators.Cubic()
   
   







