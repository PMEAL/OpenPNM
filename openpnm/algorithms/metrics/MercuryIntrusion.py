import numpy as np
from openpnm.utils import logging
from openpnm.phases import Mercury
from openpnm.physics import GenericPhysics
from openpnm.algorithms import Porosimetry, GenericAlgorithm
from openpnm import models
from openpnm import topotools
logger = logging.getLogger(__name__)


class MercuryIntrusion(GenericAlgorithm):
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

    def __init__(self, network=None, project=None, settings={}, name=None,
                 **kwargs):
        if project is None:
            project = network.project
        hg = Mercury(network=network)
        super().__init__(network=network, project=project, **kwargs)
        op = Porosimetry(project=project, phase=hg)
        self.settings['mip'] = op.name
        mod = models.physics.capillary_pressure.washburn
        for geom in project.geometries().values():
            phys = GenericPhysics(network=network, phase=hg, geometry=geom)
            phys.add_model(propname='throat.entry_pressure', model=mod)
        if not project.geometries():
            hg.add_model(propname='throat.entry_pressure', model=mod)
        topotools.find_surface_pores(network=network)
        op.set_inlets(pores=network.pores('surface'))
        logger.info('Running MIP simulation')
        op.run()
        self.update(op)

    def _set_snwp_data(self, data):
        self._snwp_data = np.array(data)

    def _get_snwp_data(self):
        if hasattr(self, '_snwp_data'):
            return self._snwp_data
        else:
            logger.error('Pc data has not been provided')

    snwp_data = property(fget=_get_snwp_data, fset=_set_snwp_data)

    def _set_pc_data(self, data):
        self._pc_data = np.array(data)

    def _get_pc_data(self):
        if hasattr(self, '_pc_data'):
            return self._pc_data
        else:
            logger.error('Pc data has not been provided')

    pc_data = property(fget=_get_pc_data, fset=_set_pc_data)

    def plot_intrusion_curve(self, fig=None):
        proj = self.project
        op = proj[self.settings['mip']]
        fig = op.plot_intrusion_curve(fig=fig)
        ax = fig.gca()
        x = self.pc_data
        y = self.snwp_data
        if (x is not None) and (y is not None):
            ax.plot(x, y, 'r*-')
        return fig
