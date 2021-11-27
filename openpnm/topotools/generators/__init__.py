r"""
================================================
Generators (:mod:`openpnm.topotools.generators`)
================================================

This module contains a selection of functions that deal specifically with
generating sufficient information that can be turned into an openpnm network.

"""

from .cubic import cubic
from .delaunay import delaunay
from .gabriel import gabriel
from .voronoi import voronoi
from .voronoi_delaunay_dual import voronoi_delaunay_dual
from .template import cubic_template
from .fcc import fcc
from .bcc import bcc
