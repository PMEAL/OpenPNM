import logging
import numpy as np
from openpnm.core import Domain
from openpnm.utils import Workspace
from openpnm.utils import Docorator, SettingsAttr
from numpy import ones
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
    This class produces a blank-slate object with no pore-scale models for
    calculating any thermophysical properties.  Users must add models and
    specify parameters for all the properties they require.

    Parameters
    ----------
    %(Base.parameters)s

    """

    def __init__(self, network, settings=None, **kwargs):
        self.settings = SettingsAttr(PhaseSettings, settings)
        super().__init__(network=network, settings=self.settings, **kwargs)

        # If project has a network object, adjust pore and throat array sizes
        self['pore.'+self.name] = np.ones([network.Np, ], dtype=bool)
        self['throat.'+self.name] = np.ones([network.Nt, ], dtype=bool)

        # Set standard conditions on the fluid to get started
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __getitem__(self, key):
        # Attempt at automatic interpolation if key not found
        element, prop = key.split('.', 1)
        if key not in self.keys():
            L = ['pore', 'throat']
            L.remove(element)
            elem = L[0]
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

    @property
    def phase(self):
        """A shortcut to get a handle to the associated phase (itself)."""
        return self

    @property
    def physics(self):
        """A shortcut to query the associated physics(es)."""
        return self.project.find_physics(phase=self)

    @property
    def _subdomains(self):
        return self.project.find_physics(phase=self)
