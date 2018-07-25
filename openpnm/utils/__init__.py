r"""

**openpnm.utils**

----

This module contains 2 very important classes (Project and Workspace) as well
as a number of helper classes.

"""

import logging as logging
from .misc import PrintableDict
from .misc import PrintableList
from .misc import NestedDict
from .misc import SettingsDict
from .misc import HealthDict
from .misc import flat_list
from .misc import sanitize_dict
from .misc import unique_list
from .misc import tic, toc
from .Workspace import Workspace
from .Project import Project


log_format = \
    '|%(levelname)-s|%(message)s'
logging.basicConfig(level=logging.WARNING, format=log_format)
del log_format
