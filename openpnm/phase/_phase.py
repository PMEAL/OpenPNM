import logging
import numpy as np
from openpnm.core import Domain
from openpnm.utils import Workspace, Docorator


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
    auto_interpolate : boolean
        If ``True`` then calls to a missing 'throat.<prop>' will automatically
        interpolate 'pore.<prop>', if present, and vice versa.  If ``False``
        a normal ``KeyError`` is raised.
    """
    auto_interpolate = True


@docstr.get_sections(base='Phase', sections=['Parameters'])
@docstr.dedent
class Phase(Domain):
    r"""
    This class produces an empty object with no pore-scale models for
    calculating any thermophysical properties.  Users must add models and
    specify parameters for all the properties they require.

    Parameters
    ----------
    %(Base.parameters)s

    """

    def __init__(self, network, name='phase_?', **kwargs):
        super().__init__(network=network, name=name, **kwargs)
        self.settings._update(PhaseSettings())
        self['pore.all'] = np.ones([network.Np, ], dtype=bool)
        self['throat.all'] = np.ones([network.Nt, ], dtype=bool)
        # Set standard conditions on the phase
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __getitem__(self, key):
        try:
            return super().__getitem__(key)
        except KeyError:
            pass

        try:
            return self.network[key]
        except KeyError:
            pass

        # Parse the key
        element, prop, domain = key.split('.', 1) + ['all']
        if '@' in prop:
            prop, domain = prop.split('@')

        # Start by getting locs
        if domain == 'all':
            locs = np.ones(self._count(element), dtype=bool)
        elif element + '.' + domain in self.keys():
            locs = super().__getitem__(element + '.' + domain)
        elif element + '.' + domain in self.network.keys():
            locs = self.network[element + '.' + domain]
        else:
            raise KeyError(element + '.' + domain)

        # Next get the data arrays
        if element + '.' + prop in self.keys():
            vals = super().__getitem__(element + '.' + prop)
        elif element + '.' + prop in self.network.keys():
            vals = self.network[element + '.' + prop]
        else:
            if self.settings['auto_interpolate']:
                if element == 'param':
                    raise KeyError(key)
                elif (element == 'pore') and ('throat.'+prop not in self.keys()):
                    msg = f"'throat.{prop}' not found, cannot interpolate '{element+'.'+prop}'"
                    raise KeyError(msg)
                elif (element == 'throat') and ('pore.'+prop not in self.keys()):
                    msg = f"'pore.{prop}', cannot interpolate '{element+'.'+prop}'"
                    raise KeyError(msg)
                vals = self.interpolate_data(element + '.' + prop)
            else:
                raise KeyError(key)
        return vals[locs]













