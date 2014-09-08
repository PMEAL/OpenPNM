r"""
###############################################################################
:mod:`OpenPNM.Physics` -- Pore Scale Physics Models
###############################################################################

Contents
--------
GenericPhysics: This base class is essentially an init statement that ensures
the Physics object registers itself with Phase and Network objects correctly.

Subclasses: OpenPNM includes one subclass called Standard which invokes the
typical pore scale physics models such as the Hagen-Poiseiulle model for 
hydraulic conductance through a tube.  Customized Physics classes can be 
created by adding a file to the Physics directory, which will be imported
automatically.  

Classes
-------
.. autoclass:: GenericPhysics
   :members:

.. autoclass:: Standard
   :members:

"""

#Import every file in the directory
import os as _os
dirname = _os.path.dirname(_os.path.abspath(__file__))
for item in _os.listdir(dirname):
    if item.split('.')[-1] == 'py':
        if item in ['__init__.py','__pycache__']:
            pass
        elif item[0:2] == '__':
            try:
                exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
            except:
                print('File name '+item+' does not match class name, cannot import automatically')
        else:
            exec('from . import ' + format(item.split('.')[0]))
