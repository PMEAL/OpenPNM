r"""
Miscellaneous
=============

This submodule contains models for calculating general or generic values
that are useful for all other pore-scale models, such as scaling values, or
generating an array of random numbers.

"""

# The following bits are to initialize some boilerplate docstrings for docrep
from openpnm.utils import Docorator as _doc
_docstr = _doc()
_docstr.params['models.misc.seeds'] = \
    r"""seeds : str, optional
            Name of the dictionary key on ``target`` where the array containing
            random seed values (between 0 and 1) is stored. These values are
            used to do a reverse lookup on the cdf of given distribution using
            the ``ppf`` method on the scipy.stats object.  Truncating the
            range (e.g. from 0.1 to 0.9) is often useful to prevent the
            occurance of extreme values from the long tails."""


from ._statistical_distributions import *
from ._simple_equations import *
from ._basic_math import *
from ._neighbor_lookups import *
