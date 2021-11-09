from openpnm.core import Base
from openpnm.utils import logging, Docorator, prettify_logger_message
import numpy as np
from openpnm.topotools import iscoplanar, is_fully_connected, dimensionality
from scipy.spatial import ConvexHull
from scipy.spatial import cKDTree
logger = logging.getLogger(__name__)


class GenericMetric(Base):
    r"""

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings['prefix'] = 'metric'
        super().__init__(project=project, **kwargs)
        project = self.project
        self['pore.all'] = np.ones((project.network.Np, ), dtype=bool)
        self['throat.all'] = np.ones((project.network.Nt, ), dtype=bool)
    
    
    def _calc_eff_prop(self, inlets=None, outlets=None,
                       domain_area=None, domain_length=None, rates=None,
                       prop_diff=None):
        r"""
        Calculates the effective transport through the network.
        Parameters
        ----------
        inlets : array_like
            The pores where the inlet boundary conditions were applied. If
            not given an attempt is made to infer them from the algorithm.
        outlets : array_like
            The pores where the outlet boundary conditions were applied.
            If not given an attempt is made to infer them from the
            algorithm.
        domain_area : float
            The area of the inlet and/or outlet face (which shold match)
        domain_length : float
            The length of the domain between the inlet and outlet faces
        Returns
        -------
        The effective transport property through the network
        """
        network = self.project.network
        if rates is None:
            raise Exception('You need to run the algorithm first.')
        Dx = prop_diff
        if domain_area is None:
            domain_area = self._get_domain_area(inlets=inlets, outlets=outlets)
        if domain_length is None:
            domain_length = self._get_domain_length(inlets=inlets, outlets=outlets)
        D = np.sum(rates) * domain_length / domain_area / Dx
        return np.atleast_1d(D)

    def _get_domain_area(self, inlets=None, outlets=None):
        r"""
        Determines the cross sectional area relative to the inlets/outlets.
        """
        logger.warning('Attempting to estimate inlet area...will be low')
        if dimensionality(self.network).sum() != 3:
            raise Exception('The network is not 3D, specify area manually')
        inlets = self.network.coords[inlets]
        outlets = self.network.coords[outlets]
        if not iscoplanar(inlets):
            logger.error('Detected inlet pores are not coplanar')
        if not iscoplanar(outlets):
            logger.error('Detected outlet pores are not coplanar')
        Nin = np.ptp(inlets, axis=0) > 0
        if Nin.all():
            logger.warning('Detected inlets are not oriented along a principle axis')
        Nout = np.ptp(outlets, axis=0) > 0
        if Nout.all():
            logger.warning('Detected outlets are not oriented along a principle axis')
        hull_in = ConvexHull(points=inlets[:, Nin])
        hull_out = ConvexHull(points=outlets[:, Nout])
        if hull_in.volume != hull_out.volume:
            logger.error('Inlet and outlet faces are different area')
        area = hull_in.volume  # In 2D: volume=area, area=perimeter
        return area

    def _get_domain_length(self, inlets=None, outlets=None):
        r"""
        Determines the domain length relative to the inlets/outlets.
        """
        msg = ('Attempting to estimate domain length...could be low if'
               ' boundary pores were not added')
        logger.warning(prettify_logger_message(msg))
        inlets = self.network.coords[inlets]
        outlets = self.network.coords[outlets]
        if not iscoplanar(inlets):
            logger.error('Detected inlet pores are not coplanar')
        if not iscoplanar(outlets):
            logger.error('Detected inlet pores are not coplanar')
        tree = cKDTree(data=inlets)
        Ls = np.unique(np.float64(tree.query(x=outlets)[0]))
        if not np.allclose(Ls, Ls[0]):
            logger.error('A unique value of length could not be found')
        length = Ls[0]
        return length
