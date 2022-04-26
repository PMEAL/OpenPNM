import logging
import numpy as np
from openpnm.core import Domain
from openpnm.utils import Workspace
from openpnm.utils import Docorator, SettingsAttr
import openpnm.models as mods


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()


@docstr.get_sections(base='PhaseSettings', sections=['Parameters'])
@docstr.dedent
class PhaseSettings:
    r"""
    Parameters
    ----------
    %(BaseSettings.parameters)s
    """
    prefix = 'phase'


@docstr.get_sections(base='GenericPhase', sections=['Parameters'])
@docstr.dedent
class GenericPhase(Domain):
    r"""
    This class produces an empty object with no pore-scale models for
    calculating any thermophysical properties.  Users must add models and
    specify parameters for all the properties they require.

    Parameters
    ----------
    %(Base.parameters)s

    """

    def __init__(self, network, settings=None, **kwargs):
        self.settings = SettingsAttr(PhaseSettings, settings)
        super().__init__(network=network, settings=self.settings, **kwargs)

        self['pore.all'] = np.ones([network.Np, ], dtype=bool)
        self['throat.all'] = np.ones([network.Nt, ], dtype=bool)

        # Set standard conditions on the fluid to get started
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __getitem__(self, key):
        # Attempt at automatic interpolation if key not found.  This behavior
        # could arguably be turned off/on as a setting
        element, prop = key.split('.', 1)
        if key not in self.keys():
            elem = list({'pore', 'throat'}.difference({element}))[0]
            if (elem + '.' + prop) in self.keys():
                mod = {'pore': mods.misc.from_neighbor_throats,
                       'throat': mods.misc.from_neighbor_pores}
                self.add_model(propname=key,
                               model=mod[element],
                               prop=elem + '.' + prop,
                               mode='mean')
                self.run_model(key)
        vals = super().__getitem__(key)
        return vals
