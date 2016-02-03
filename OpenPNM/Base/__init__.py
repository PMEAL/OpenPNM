r"""
###############################################################################
:mod:`OpenPNM.Base` -- Module Containing Abstract Base Class, Core Data Class,
and other backend classes used in OpenPNM
###############################################################################

.. autoclass:: OpenPNM.Base.Workspace
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: OpenPNM.Base.Core
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: OpenPNM.Base.Tools
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: OpenPNM.Base.ModelsDict
   :members:
   :undoc-members:
   :show-inheritance:

"""
import logging as logging
from .__Workspace__ import Workspace
from .__Controller__ import Controller
from .__ModelsDict__ import ModelsDict
from . import __Tools__ as Tools
from .__Core__ import Core


# Set up logging to file - see previous section for more details
log_format = \
    '%(asctime)s | %(levelname)-8s | %(name)s.%(funcName)s | %(message)s'
logging.basicConfig(level=logging.WARNING, format=log_format)
