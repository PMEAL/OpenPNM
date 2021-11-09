import numpy as np
from openpnm.utils import logging, Project, Workspace, PrintableDict
from openpnm.phases import GenericPhase
from openpnm.physics import GenericPhysics
from openpnm.algorithms import StokesFlow
from openpnm.metrics import GenericMetric
from openpnm import models
from openpnm import topotools
logger = logging.getLogger(__name__)
ws = Workspace()

default_settings = {
    'inlet': None,
    'outlet': None,
    'area': None,
    'length': None,
}

class AbsolutePermeability(GenericMetric):
    r"""
    This class works by applying 'value' boundary conditions across the
    domain to find molar flow, then using Fick's law to back-calculate
    the effective permeability of the domain.

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings.update(default_settings)
        self.settings.update(settings)
        if network is None:
            network = project.network
        if project is None:
            project = network.project
        super().__init__(network=network, project=project, **kwargs)

    def run(self):
        r"""
        Execute the diffusion simulations in the principle directions.

        """
        phase = GenericPhase(network=self.network)
        phase['pore.viscosity'] = 1.0
        phase['throat.viscosity'] = 1.0
        mod = models.physics.hydraulic_conductance.hagen_poiseuille
        for geom in self.project.geometries().values():
            phys = GenericPhysics(network=self.network,
                                  phase=phase, geometry=geom)
            phys.add_model(propname='throat.hydraulic_conductance', model=mod)
        inlet = self.network.pores(self.settings['inlet'])
        outlet = self.network.pores(self.settings['outlet'])
        perm = StokesFlow(network=self.project.network, phase=phase)
        perm.set_value_BC(pores=inlet, values= 1)
        perm.set_value_BC(pores=outlet, values=0)
        perm.run()
        phase.update(perm.results())
        K= self._calc_eff_prop(inlets=inlet, outlets=outlet,
                       domain_area=None, domain_length=None, rates=perm.rate(pores=inlet),
                       prop_diff=1)
        return K
