r"""
Generators
==========

This module contains a selection of functions that generate various network
topologies

"""

from ._cubic import cubic
from ._delaunay import delaunay
from ._gabriel import gabriel
from ._voronoi import voronoi
from ._voronoi_delaunay_dual import voronoi_delaunay_dual
from ._template import cubic_template
from ._fcc import fcc
from ._bcc import bcc
from . import tools
