r"""
*******************************************************************************
:mod:`OpenPNM.Geometry` -- Classes related to the creation of pore and throat geometry
*******************************************************************************

.. module:: OpenPNM.Geometry

Contents
--------
Contains methods for applying pore and throat geometry

Classes
-------
    
.. autoclass:: GenericGeometry
   :members:
   :undoc-members:
   :show-inheritance:
   
.. autoclass:: Stick_and_Ball
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: Boundary
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.pore_seed
   :members:
   :undoc-members:
   :show-inheritance:
   
.. automodule:: OpenPNM.Geometry.pore_diameter
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.pore_volume
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.throat_seed
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.throat_diameter
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.throat_volume
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.throat_length
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.throat_area
   :members:
   :undoc-members:
   :show-inheritance:

.. automodule:: OpenPNM.Geometry.throat_surface_area
   :members:
   :undoc-members:
   :show-inheritance:
   
"""
#Import every file in the directory, including both classes and methods
#import os,sys
#dir = os.path.dirname(os.path.abspath(__file__))
#for item in os.listdir(dir):
#    if item.split('.')[-1] == 'py':
#        if item == '__init__.py':
#            pass
#        elif item[0:2] == '__':
#            exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
#        else:
#            exec('from . import ' + format(item.split('.')[0]))
        
from . import pore_diameter
from . import pore_seed
from . import pore_volume
from . import pore_centroid
from . import throat_centroid
from . import throat_diameter
from . import throat_length
from . import throat_seed
from . import throat_volume
from . import throat_vector
from . import throat_area
from . import throat_surface_area

