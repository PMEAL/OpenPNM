r"""
Collection of pore-scale models
===============================

"""


# The following bits are to initialize some boilerplate docstrings
from openpnm.utils import Docorator as _doc
from matplotlib.docstring import Substitution as _sub


_docstr = _doc()
_docstr.params['models.target.parameters'] = \
    r"""target : OpenPNM Base object
                Object with which this model is associated. This controls
                the length of the calculated array, and also provides access
                to other necessary properties."""


_doctxt = _sub(
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


from . import misc
from . import network
from . import geometry
from . import phase
from . import physics
from . import collections
