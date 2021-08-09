r"""

**openpnm.utils**

----

This module contains 2 very important classes (Project and Workspace) as well
as a number of helper classes.

"""

import logging as logging
from .misc import Docorator
from .misc import PrintableDict
from .misc import PrintableList
from .misc import NestedDict
from .misc import SubDict
from .misc import SettingsDict
from .misc import GenericSettings
from .misc import HealthDict
from .misc import flat_list
from .misc import sanitize_dict
from .misc import unique_list
from .misc import tic, toc
from .misc import is_symmetric
from .misc import nbr_to_str
from .misc import conduit_dict_to_array
from .misc import conduit_array_to_dict
from .misc import prettify_logger_message
from .Workspace import Workspace
from .Project import Project


# You can add info to the logger message by inserting the desired %(item)
# For a list of available items see:
# https://docs.python.org/3/library/logging.html#logrecord-attributes
# NOTE: If the calling locations appears as 'root' it's because the logger
# was not given a name in a file somewhere.  A good option is __name__.
log_format = \
'-' * 60 + '\n\
%(levelname)-11s: %(message)s \n\
SOURCE     : %(name)s.%(funcName)s \n\
TIME STAMP : %(asctime)s\
\n' + '-' * 60

logging.basicConfig(level=logging.WARNING, format=log_format)
del log_format
