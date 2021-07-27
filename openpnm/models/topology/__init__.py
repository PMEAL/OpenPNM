r"""

**openpnm.models.network**

----

This submodule contains models for calculating topological properties of
networks

"""

from .topology import coordination_number
from .topology import total_length
from .health import isolated_pores
from .health import oversize_throats
from .health import bidirectional_throats
from .health import headless_throats
