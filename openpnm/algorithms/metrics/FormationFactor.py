import numpy as np
from openpnm.utils import logging, Project, Workspace, PrintableDict
from openpnm.phases import GenericPhase
from openpnm.physics import GenericPhysics
from openpnm.algorithms import FickianDiffusion
from openpnm.algorithms.metrics import GenericMetric
from openpnm import models
from openpnm import topotools
logger = logging.getLogger(__name__)
ws = Workspace()


class FormationFactor(GenericMetric):
    r"""
    Nothing yet

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-5)
    >>> geo = op.geometry.StickAndBall(network=pn)

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings.update({'inlets': {'x': 'left',
                                         'y': 'front',
                                         'z': 'top'},
                              'outlets': {'x': 'right',
                                          'y': 'back',
                                          'z': 'bottom'},
                              'areas': {'x': None,
                                        'y': None,
                                        'z': None},
                              'lengths': {'x': None,
                                          'y': None,
                                          'z': None}})

        self.results = PrintableDict()
        self.results._value = "Formation Factor"
        self.results._key = "Direction"
        if network is None:
            network = project.network
        if project is None:
            project = network.project
#        project = ws.copy_project(network.project)
#        keep = [network] + list(project.geometries().values())
#        for i in project:
#            if i not in keep:
#                project.purge_object(i)
        super().__init__(network=network, project=project, **kwargs)

    def run(self):
        r"""
        Execute the diffusion simulations in the principle directions

        Notes
        -----
        The default settings
        """
        phase = GenericPhase(network=self.network)
        phase['pore.diffusivity'] = 1.0
        phase['throat.diffusivity'] = 1.0
        mod = models.physics.diffusive_conductance.ordinary_diffusion
        for geom in self.project.geometries().values():
            phys = GenericPhysics(network=self.network,
                                  phase=phase, geometry=geom)
            phys.add_model(propname='throat.diffusive_conductance', model=mod)
        for bcs in self.settings['inlets'].keys():
            Diff = FickianDiffusion(network=self.project.network, phase=phase)
            Pin = self.network.pores(self.settings['inlets'][bcs])
            Diff.set_value_BC(pores=Pin, values=1.0)
            Pout = self.network.pores(self.settings['outlets'][bcs])
            Diff.set_value_BC(pores=Pout, values=0.0)
            Diff.run()
            A = self.settings['areas'][bcs]
            if A is None:
                A = Diff._get_domain_area()
            L = self.settings['lengths'][bcs]
            if L is None:
                L = Diff._get_domain_length()
            R = Diff.rate(pores=Pin)
            Deff = R*L/A  # Conc gradient and diffusivity were both unity
            self.results[bcs] = 1/Deff[0]
        print(self.results)

    def set_inlets(self, direction, label):
        self.settings['inlets'][direction] = label

    def set_outlets(self, direction, label):
        self.settings['outlets'][direction] = label

    def set_area(self, direction, area):
        self.settings['areas'][direction] = area

    def set_length(self, direction, length):
        self.settings['lengths'][direction] = length
