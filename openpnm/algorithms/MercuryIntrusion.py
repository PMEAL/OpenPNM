from openpnm.utils import logging
from openpnm.phases import Mercury
from openpnm.physics import GenericPhysics
from openpnm.algorithms import Porosimetry
from openpnm import models
from openpnm import topotools
logger = logging.getLogger(__name__)


class MercuryIntrusion(Porosimetry):
    r"""
    A ready-make Mercury Intrusion Porosimetry algorithm

    This class accepts a pre-existing Project, including a Network and any
    associated Geometries, then internally creates all the necessary objects
    for running an MIP simulation.  Specifically this means creating an Hg
    Phase object, and adding the Washburn capillary pressure model.
    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        if project is None:
            project = network.project
        hg = Mercury(network=network)
        kwargs['phase'] = hg
        super().__init__(project=project, **kwargs)
        mod = models.physics.capillary_pressure.washburn
        for geom in project.geometries().vaues():
            phys = GenericPhysics(network=network, phase=hg, geometry=geom)
            phys.add_model(propname='throat.entry_pressure', model=mod)
        self.setup(phase=hg)
        topotools.find_surface_pores(network=network)
        self.set_inlets(pores=network.pores('surface'))
