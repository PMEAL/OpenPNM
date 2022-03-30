r"""
Geometry
--------

This submodule contains pore-scale models that calculate geometrical
properties. These models are to be added to a Geometry object.

"""

# The following bits are to initialize some boilerplate docstrings for docrep
from openpnm.utils import Docorator as _doc
_docstr = _doc()
_docstr.params['models.geometry.pdia'] = \
    r"""pore_diameter : str
            Name of the dictionary key on ``target`` where the array containing
            pore diameter values is stored"""
_docstr.params['models.geometry.tdia'] = \
    r"""throat_diameter : str
            Name of the dictionary key on ``target`` where the array containing
            pore diameter values is stored"""
_docstr.params['models.geometry.pvol'] = \
    r"""pore_volume : str
            Name of the dictionary key on ``target`` where the array containing
            pore volume values is stored"""
_docstr.params['models.geometry.tvol'] = \
    r"""throat_volume : str
            Name of the dictionary key on ``target`` where the array containing
            throat volume values is stored"""
_docstr.params['models.geometry.tlen'] = \
    r"""throat_length : str
            Name of the dictionary key on ``target`` where the array containing
            throat length values is stored"""
_docstr.params['models.geometry.tarea'] = \
    r"""throat_area : str
            Name of the dictionary key on ``target`` where the array containing
            throat area values is stored"""

# Perform imports
from . import pore_size
from . import pore_seed
from . import pore_volume
from . import pore_surface_area
from . import pore_cross_sectional_area
from . import throat_cross_sectional_area
from . import throat_seed
from . import throat_size
from . import throat_length
from . import throat_perimeter
from . import throat_surface_area
from . import throat_volume
from . import throat_capillary_shape_factor
from . import throat_centroid
from . import throat_vector
from . import hydraulic_size_factors
from . import diffusive_size_factors
from . import conduit_lengths
