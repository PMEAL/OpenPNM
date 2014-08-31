r"""
###############################################################################
:mod:`OpenPNM.Phases` -- Phase Property Estimation Methods
###############################################################################

.. module:: OpenPNM.Phases

Contents
--------
This module contains subclasses for producing specific phases


Classes
-------

.. autoclass:: GenericPhase
   :members:

.. autoclass:: Air
   :members:

.. autoclass:: Water
   :members:
   
.. autoclass:: Mercury
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
            exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
        else:
            exec('from . import ' + format(item.split('.')[0]))
