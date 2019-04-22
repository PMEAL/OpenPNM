from openpnm.core import Base
from openpnm.utils import logging
logger = logging.getLogger(__name__)

# Set some default settings
def_set = {'pore_volume': 'pore.volume',
           'throat_volume': 'throat.volume'}


class Porosity(Base):
    r"""
    This class implements steady-state linear transport calculations

    Parameters
    ----------
    network : OpenPNM Network object
        The Network with which this algorithm is associated

    project : OpenPNM Project object, optional
        A Project can be specified instead of ``network``

    """

    def __init__(self, project=None, network=None, settings={},
                 **kwargs):
        # Apply default settings
        self.settings.update(def_set)
        # Overwrite any given in init
        self.settings.update(settings)
        # If network given, get project, otherwise let parent class create it
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)

    def setup(self, pore_volume=None, throat_volume=None):
        if pore_volume is not None:
            self.settings['pore_volume'] = 'pore.volume'
        if throat_volume is not None:
            self.settings['throat_volume'] = 'throat.volume'

    def find_intersection_volume(self):
        return 0

    def porosity(self):
        Vvoid = self.get_pore_volume() + self.get_throat_volume() \
                - self.find_intersection_volume()
        phi = Vvoid/self.bulk_volume
        return phi

    def get_pore_volume(self):
        Pvol = self.project.network[self.settings['pore_volume']].sum()
        return Pvol

    def get_throat_volume(self):
        Tvol = self.project.network[self.settings['throat_volume']].sum()
        return Tvol

    def _get_Vbulk(self):
        if not hasattr(self, '_Vbulk'):
            self._Vbulk = None
        return self._Vbulk

    def _set_Vbulk(self, val):
        self._Vbulk = val

    bulk_volume = property(fget=_get_Vbulk, fset=_set_Vbulk)
