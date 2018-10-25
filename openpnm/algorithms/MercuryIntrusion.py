from openpnm.utils import logging
from openpnm.phases import Mercury
from openpnm.physics import GenericPhysics
from openpnm.algorithms import Porosimetry
from openpnm import models
from openpnm import topotools
logger = logging.getLogger(__name__)


class MercuryIntrusion(Porosimetry):
    r"""
    A ready-made Mercury Intrusion Porosimetry algorithm

    This class accepts a pre-existing Project, including a Network and any
    associated Geometries, then internally creates a new Mercury Phase object
    and a new Physics for each Geometry associated with the Network.  Finally,
    the Washburn capillary pressure model is added to each Physics object (
    or to the Mercury object if there are no Geometries defined).

    Notes
    -----
    The simulation is automatically run with all faces treated as inlets.  The
    results can be plotted using `plot_intrusion_data`, and numerical data
    can be obtained with `get_intrusion_data`.
    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        if project is None:
            project = network.project
        hg = Mercury(network=network)
        kwargs['phase'] = hg
        super().__init__(project=project, **kwargs)
        mod = models.physics.capillary_pressure.washburn
        for geom in project.geometries().values():
            phys = GenericPhysics(network=network, phase=hg, geometry=geom)
            phys.add_model(propname='throat.entry_pressure', model=mod)
        if not project.geometries():
            hg.add_model(propname='throat.entry_pressure', model=mod)
        self.setup(phase=hg)
        topotools.find_surface_pores(network=network)
        self.set_inlets(pores=network.pores('surface'))
        logger.info('Running MIP simulation')
        self.run()
