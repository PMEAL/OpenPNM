r"""
###############################################################################
:mod:`OpenPNM.Network`: Classes related the creation of network topology
###############################################################################

.. module:: OpenPNM.Network

Contents
--------
GenericNetwork: Contains many classes for working with the topology of the
networks

Subclasses: Inherit from GenericNetwork, and contain additional methods for
actually generating topology.  Users can create thrir own topology generators
by adding a file to the Network directory.  

Classes
-------

.. autoclass:: GenericNetwork
   :members:

.. autoclass:: Cubic
   :members:

.. autoclass:: Delaunay
   :members:


"""
import OpenPNM.Base
#Import every file in the directory
import os as _os
dir = _os.path.dirname(_os.path.abspath(__file__))
for item in sorted(_os.listdir(dir)):
    if item.endswith('.py'):
        if item == '__init__.py':
            pass
        elif item.startswith('__'):
            exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
        else:
            exec('from . import ' + format(item.split('.')[0]))
            
#Manually added imports
