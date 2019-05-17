r"""

**openpnm.models.physics**

----

This submodule contains models for calculating properties related to physical
transport processes, including conductances, reaction rates, and capillary
effects

"""

from . import capillary_pressure
from . import diffusive_conductance
from . import electrical_conductance
from . import thermal_conductance
from . import hydraulic_conductance
from . import dispersive_conductance
from . import multiphase
from . import generic_source_term
from . import flow_shape_factors
from . import poisson_shape_factors
from . import meniscus
from . import ionic_conductance
from . import ad_dif_mig_conductance
from . import ad_dif_conductance
