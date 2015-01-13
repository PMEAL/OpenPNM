r"""
###############################################################################
:mod:`OpenPNM.Base` -- Abstract Base Class, and Core Data Class
###############################################################################

.. autoclass:: OpenPNM.Base.Controller
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: OpenPNM.Base.Base
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: OpenPNM.Base.Core
   :members:
   :undoc-members:
   :show-inheritance:

"""

import logging as logging
# set up logging to file - see previous section for more details
logging.basicConfig(level=logging.WARNING,
                    format='%(asctime)s | %(levelname)-8s | %(name)s.%(funcName)s | %(message)s',
                    )

from .__Controller__ import Controller
from .__Core__ import Core
