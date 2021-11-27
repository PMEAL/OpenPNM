r"""
Models
======

"""


# %% The following bits are to initialize some boilerplate docstrings for docrep
from openpnm.utils import Docorator as _doc
_docstr = _doc()
_docstr.params['models.target.parameters'] = \
    r"""target : OpenPNM Base object
            Object with which this model is associated. This controls
            the length of the calculated array, and also provides access to
            other necessary properties."""


# %%


__all__ = [
    'misc',
    'network',
    'geometry',
    'phases',
    'physics',
    ]

from . import misc
from . import network
from . import geometry
from . import phases
from . import physics
