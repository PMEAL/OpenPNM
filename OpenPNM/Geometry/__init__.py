r"""
###############################################################################
:mod:`OpenPNM.Geometry` -- Classes related to the creation of pore and throat geometry
###############################################################################

.. module:: OpenPNM.Geometry

Contents
--------
GenericGeometry: The init statement performs the important tasks of registering
itself with the network and creating data dictionaries of the necessary size

Subclasses: These can be custom made by users to represent specfic geometries.
It is necessary for their init's to call the init of the GenericGeometry class
in order to be properly instantiated.  OpenPNM comes with a few basic 
pre-written materials.  New classes added to this directory will be 
automatically imported and available.

Classes
-------

.. autoclass:: GenericGeometry
   :members:

.. autoclass:: Stick_and_Ball
   :members:

.. autoclass:: Cube_and_Cuboid
   :members:


   
"""

#Import every file in the directory
import os as _os
dir = _os.path.dirname(_os.path.abspath(__file__))
for item in _os.listdir(dir):
    if item.split('.')[-1] == 'py':
        if item == '__init__.py':
            pass
        elif item[0:2] == '__':
            try:
                exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
            except:
                print('File name '+item+' does not match class name, cannot import automatically')
        else:
            exec('from . import ' + format(item.split('.')[0]))
        
