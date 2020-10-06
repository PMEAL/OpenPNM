r"""

**openpnm.models.misc**

----

This submodule contains models for calculating general or generic values
that are useful for all other pore-scale models, such as scaling values, or
generating an array of random numbers.

"""

from .statistical_distributions import generic_distribution, normal, weibull, random
from .simple_equations import linear, polynomial, generic_function
from .basic_math import constant, product, scaled, clip, normalize, summation
from .basic_math import fraction, invert, blank
from .neighbor_lookups import from_neighbor_throats, from_neighbor_pores
