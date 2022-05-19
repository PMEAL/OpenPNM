import logging
import numpy as np
from openpnm.core import Domain
from openpnm.utils import Workspace
from openpnm.utils import Docorator, SettingsAttr


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
        try:
            # return super().__getitem__(key)
            try:
                return super().__getitem__(key)
            except KeyError:
                return self.network[key]
        except KeyError:
            # But also need to handle the @ look-ups somehow
            if '@' in key:
                raise Exception('Interpolation of @domain values is not supported')
            # Before interpolating, ensure other prop is present, to avoid
            # infinite recurrsion
            element, prop = key.split('.', 1)
            if (element == 'pore') and ('throat.'+prop not in self.keys()):
                raise KeyError(f"Cannot interpolate '{element+'.'+prop}' without 'throat.{prop}'")
            elif (element == 'throat') and ('pore.'+prop not in self.keys()):
                raise KeyError(f"Cannot interpolate '{element+'.'+prop}' without 'pore.{prop}'")
            vals = self.interpolate_data(element + '.' + prop)
            return vals
