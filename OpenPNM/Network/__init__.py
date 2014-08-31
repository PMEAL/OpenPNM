r"""
###############################################################################
:mod:`OpenPNM.Network`: Classes related the creation of network topology
###############################################################################

.. module:: OpenPNM.Network

Contents
--------
GenericNetwork: Contains many classes for working with the topology of the
networks

Subclasses: Inherit fro GenericNetwork, and contain additional methods for
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
            
#Manually added imports
from .__Import__ import MatFile