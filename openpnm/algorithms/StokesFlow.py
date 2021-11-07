import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator
from traits.api import Str
logger = logging.getLogger(__name__)
docstr = Docorator()


class StokesFlowSettings:
    quantity = 'pore.pressure'
    conductance = 'throat.hydraulic_conductance'


class StokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings._update(StokesFlowSettings, docs=False)
        self.settings._update(settings)

    def calc_effective_permeability(self, inlets=None, outlets=None,
                                    domain_area=None, domain_length=None):
        r"""
        This calculates the effective permeability in this linear transport
        algorithm.

        Parameters
        ----------
        inlets : array_like
            The pores where the inlet pressure boundary conditions were
            applied.  If not given an attempt is made to infer them from the
            algorithm.

        outlets : array_like
            The pores where the outlet pressure boundary conditions were
            applied.  If not given an attempt is made to infer them from the
            algorithm.

        domain_area : scalar, optional
            The area of the inlet (and outlet) boundary faces.  If not given
            then an attempt is made to estimate it, but it is usually
            underestimated.

        domain_length : scalar, optional
            The length of the domain between the inlet and outlet boundary
            faces.  If not given then an attempt is made to estimate it, but it
            is usually underestimated.

        Notes
        -----
        The area and length of the domain are found using the bounding box
        around the inlet and outlet pores which do not necessarily lie on the
        edge of the domain, resulting in underestimation of sizes.
        """
        phase = self.project.phases()[self.settings['phase']]
        d_normal = self._calc_eff_prop(inlets=inlets, outlets=outlets,
                                       domain_area=domain_area,
                                       domain_length=domain_length)
        K = d_normal * np.mean(phase['pore.viscosity'])
        return K
