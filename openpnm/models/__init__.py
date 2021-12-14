r"""
Collection of pore-scale models for manipulating data
=====================================================

"""


# %% The following bits are to initialize some boilerplate docstrings
from openpnm.utils import Docorator
from matplotlib.docstring import Substitution


_docstr = Docorator()
_docstr.params['models.target.parameters'] = \
    r"""target : OpenPNM Base object
                Object with which this model is associated. This controls
                the length of the calculated array, and also provides access
                to other necessary properties."""


_doctxt = Substitution(
    dict_blurb=\
    r"""Name of the dictionary key on ``target`` pointing to the
    array containing values of """,
    target_blurb=\
    r"""target : OpenPNM Base object
            Object with which this model is associated. This controls
            the length of the calculated array, and also provides access
            to other necessary properties.""",
    return_arr=\
    r"""values : ndarray
            A numpy ndarray containing the computed values of """,
)


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
