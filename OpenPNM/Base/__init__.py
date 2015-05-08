r"""
###############################################################################
:mod:`OpenPNM.Base` -- Abstract Base Class, and Core Data Class
###############################################################################

.. autoclass:: OpenPNM.Base.Controller
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

"""
import logging as logging
from .__Controller__ import Controller
from .__ModelsDict__ import ModelsDict
from . import __Tools__ as Tools
from .__Core__ import Core

# Set up logging to file - see previous section for more details
log_format = \
    '%(asctime)s | %(levelname)-8s | %(name)s.%(funcName)s | %(message)s'
logging.basicConfig(level=logging.WARNING, format=log_format)
