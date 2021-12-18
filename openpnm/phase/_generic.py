from openpnm.core import Base, ModelsMixin, ParamMixin, LabelMixin
from openpnm.utils import Workspace, logging
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
class GenericPhase(ParamMixin, Base, ModelsMixin, LabelMixin):
    r"""
    This generic class is meant as a starter for custom Phase objects

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
        self['pore.all'] = ones((network.Np, ), dtype=bool)
        self['throat.all'] = ones((network.Nt, ), dtype=bool)

        # Set standard conditions on the fluid to get started
        self['pore.temperature'] = 298.0
        self['pore.pressure'] = 101325.0

    def __getitem__(self, key):
        # If the key is a just a numerical value, the kick it directly back
        # This allows one to do either value='pore.blah' or value=1.0
        if isinstance(key, (int, float, bool, complex)):
            return key
        element, prop = key.split('.', 1)
        # Deal with special keys first
        if prop == '_id':
            net = self.project.network
            return net[f"{element}._id"]
        if prop == self.name:
            return self[f"{element}.all"]
        # An attempt at automatic interpolation if key not found
        if key not in self.keys():
            not_el = list(set(['pore', 'throat']).difference(set([element])))[0]
            if (not_el + '.' + prop) in self.keys():
                mod = {'pore': mods.misc.from_neighbor_throats,
                       'throat': mods.misc.from_neighbor_pores}
                self.add_model(propname=key,
                               model=mod[element],
                               prop=not_el + '.' + prop,
                               mode='mean')
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
