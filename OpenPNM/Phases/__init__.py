r"""
###############################################################################
:mod:`OpenPNM.Phases` -- Phase Property Estimation Methods
###############################################################################

Contents
--------
GenericPhase: The basic class which defines how a Phase is instantiated.  It
also has a few specific method for querying the health of mixtures or physics
objects.

Subclasses: OpenPNM includes a few pre-written subclasses that describe the 
most commonly used materials, like Air, Water and Mercury.  Creating a custom
Phase subclass simply requires placing a file in the Phases directory and it 
will be automatically loaded.  

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
            try:
                exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
            except:
                print('File name '+item+' does not match class name, cannot import automatically')
        else:
            exec('from . import ' + format(item.split('.')[0]))
