r"""
*******************************************************************************
:mod:`OpenPNM.Algorithms` -- Algorithms on Networks
*******************************************************************************

.. module:: OpenPNM.Algorithms

Contents
--------
This submodule contains algorithms for performing simulations on pore networks

.. autoclass:: GenericAlgorithm
   :members:

.. autoclass:: OrdinaryPercolation
   :members:

.. autoclass:: InvasionPercolation
   :members:
   
.. autoclass:: InvasionPercolationForImbibition
   :members:

.. autoclass:: FickianDiffusion
   :members:

.. autoclass:: StokesFlow
   :members:
   
.. autoclass:: OhmicConduction
   :members:
   
.. autoclass:: FourierConduction
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
