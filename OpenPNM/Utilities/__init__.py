r"""
*******************************************************************************
:mod:`OpenPNM.Utilities` -- IO, geometry tools and other functions
*******************************************************************************

.. module:: OpenPNM.Utilities 


"""

##Import every file in the directory
#import os as _os
#dirname = _os.path.dirname(_os.path.abspath(__file__))
#for item in _os.listdir(dirname):
#    if item.split('.')[-1] == 'py':
#        if item in ['__init__.py','__pycache__']:
#            pass
#        elif item[0:2] == '__':
#            exec('from .' + item.split('.')[0] + ' import ' + item.split('__')[1])
#        else:
#            exec('from . import ' + format(item.split('.')[0]))

from . import transformations
from . import misc
from . import IO
            