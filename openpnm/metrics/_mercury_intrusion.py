import logging
from copy import deepcopy
from openpnm.phase import Mercury
from openpnm.algorithms import Drainage
from openpnm.utils import Project
from openpnm import models
from openpnm import topotools
logger = logging.getLogger(__name__)


__all__ = [
    'MercuryIntrusionPorosimetry',
]


class MercuryIntrusionPorosimetry(Project):
    r"""
    A ready-made Mercury Intrusion Porosimetry algorithm for probing the pore-
    size distribution of a network.

    This metric accepts a network with predefined geometry properties,
    then creates a new mercury phase object and applies the Washburn
    capillary pressure model. It then applies inlet boundary conditions to
    *all* surface pores (found using ``topotools.find_surface_pores``) and
    runs the algorithm automatically.

    """

    def __init__(self, network, **kwargs):
        super().__init__(network=network, **kwargs)
        net = deepcopy(network)
        self.append(net)
        hg = Mercury(network=net)
        mod = models.physics.capillary_pressure.washburn
        hg.add_model(propname='throat.entry_pressure', model=mod)
        self._drn = Drainage(network=net, phase=hg)
        topotools.find_surface_pores(network=net)
        self._drn.set_inlets(pores=net.pores('surface'))
        self._drn.run()

    @property
    def pc_data(self):
        data = self.get_intrusion_data()
        return data.pc

    @property
    def snwp_data(self):
        data = self.get_intrusion_data()
        return data.snwp

    def get_intrusion_data(self, pressures=None):
        return self._drn.pc_curve(pressure=pressures)
