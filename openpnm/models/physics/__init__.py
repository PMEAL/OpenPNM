r"""

**openpnm.models.physics**

----

This submodule contains models for calculating properties related to physical
transport processes, including conductances, reaction rates, and capillary
effects

"""

from . import ad_dif_mig_conductance
from . import ad_dif_conductance
from . import diffusive_conductance
from . import dispersive_conductance
from . import electrical_conductance
from . import hydraulic_conductance
from . import ionic_conductance
from . import thermal_conductance
from . import source_terms
from . import source_terms as generic_source_term
from . import capillary_pressure
from . import meniscus
from . import multiphase
from . import poisson_shape_factors
from . import flow_shape_factors
