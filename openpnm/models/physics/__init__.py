r"""
Physics
-------

This submodule contains models for calculating properties related to physical
transport processes, including conductances, reaction rates, and capillary
effects

"""

# %% The following bits are to define some boilerplate docstrings for docrep
from openpnm.utils import Docorator as _doc


_docstr = _doc()


@_docstr.get_sections(base='models.phys', sections=['Parameters'])
def _dummy():
    r"""
    pore_temperature : str
        Name of the dictionary key on ``target`` where the array of
        pore temperature
    pore_pressure : str
        Name of the dictionary key on ``target`` where the array containing
        pore_pressure values is stored
    pore_diffusivity : str
        Name of the dictionary key on ``target`` where the array containing
        pore diffusivity values is stored
    throat_diffusivity : str
        Name of the dictionary key on ``target`` where the array containing
        throat diffusivity values is stored
    pore_viscosity : str
        Name of the dictionary key on ``target`` where the array containing
        pore viscosity values is stored
    throat_viscosity : str
        Name of the dictionary key on ``target`` where the array containing
        throat viscosity values is stored
    """
    # This is a dummy function as a place to write boilerplate docstring
    # components for use in the models.physics module
    pass


from . import ad_dif_mig_conductance
from . import ad_dif_conductance
from . import diffusive_conductance
from . import electrical_conductance
from . import hydraulic_conductance
from . import ionic_conductance
from . import thermal_conductance
from . import source_terms
from . import source_terms as generic_source_term
from . import capillary_pressure
from . import meniscus
from . import multiphase
