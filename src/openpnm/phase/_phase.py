import logging
import numpy as np
from openpnm.core import Domain
from openpnm.utils import Workspace, Docorator


docstr = Docorator()
logger = logging.getLogger(__name__)
ws = Workspace()


__all__ = [
    'Phase',
]


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
    default_domain = 'all'


@docstr.get_sections(base='Phase', sections=['Parameters'])
@docstr.dedent
class Phase(Domain):
    r"""
    This class produces an empty object with no pore-scale models for
    calculating any thermophysical properties. Users must add models and
    specify parameters for all the properties they require.

    Parameters
    ----------
    %(Base.parameters)s

    """

    def __init__(self, network, name='phase_?', **kwargs):
        super().__init__(network=network, name=name, **kwargs)
        self.settings._update(PhaseSettings())
        # Set standard conditions on the phase
        self['pore.all'] = np.ones([network.Np, ], dtype=bool)
        self['throat.all'] = np.ones([network.Nt, ], dtype=bool)
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __setitem__(self, key, value):
        if '@' not in key:
            super().__setitem__(key, value)
        else:
            # Deal with the fact that the label might only exist on the network
            propname, domain = key.split('@')
            element, prop = propname.split('.', 1)
            try:  # Fetch array from self if present
                temp = self[element + '.' + prop]
            except KeyError:  # Otherwise create it
                temp = self._initialize_empty_array_like(value, element)
                self[element + '.' + prop] = temp
            # Insert values into masked locations
            mask = self.project._get_locations(element + '.' + domain)
            temp[mask] = value

    def __getitem__(self, key):
        try:  # If key exists, just get it
            return super().__getitem__(key)
        except KeyError:
            pass

        try:  # Allow look-up from network mostly for label/domain info
            return self.network[key]
        except (KeyError, AttributeError):
            pass

        # Parse the key
        element, prop = key.split('.', 1)
        if '@' in prop:
            prop, domain = prop.split('@')
        else:
            domain = 'all'

        # Get params directly if appropriate
        if element == 'param':
            return self.params[prop]

        # Next get the data arrays, this is the case if @ notation was used
        if element + '.' + prop in self.keys():
            vals = super().__getitem__(element + '.' + prop)
        else:  # If above are not triggered then try to interpolate
            try:
                if self.settings['auto_interpolate']:
                    if (element == 'pore') and ('throat.'+prop not in self.keys()):
                        raise KeyError(key)
                    elif (element == 'throat') and ('pore.'+prop not in self.keys()):
                        raise KeyError(key)
                    vals = self.interpolate_data(element + '.' + prop)
                else:
                    raise KeyError(key)
            except AttributeError:
                raise KeyError(key)

        # Finally get locs
        if domain == 'all':
            locs = np.ones(self._count(element), dtype=bool)
        else:
            locs = self.project._get_locations(element + '.' + domain)

        return vals[locs]
