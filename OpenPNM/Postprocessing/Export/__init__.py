r"""
*******************************************************************************
:mod:`OpenPNM.Postprocessing.Export` -- 
*******************************************************************************

.. module:: OpenPNM.Postprocessing.Export

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