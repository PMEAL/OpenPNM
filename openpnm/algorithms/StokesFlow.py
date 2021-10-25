import numpy as np
import scipy as sp
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, GenericSettings, Docorator
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='StokesFlowSettings',
                     sections=['Parameters'])
@docstr.dedent
class StokesFlowSettings(GenericSettings):
    r"""

    Parameters
    ----------
    %(GenericTransportSettings.parameters)s
    quantity : str (default = 'pore.pressure')
        The name of the physical quantity to be calculated
    conductance : str (default = 'throat.hydraulic_conductance')
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.

    Other Parameters
    ----------------

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """
    quantity = 'pore.pressure'
    conductance = 'throat.hydraulic_conductance'


class StokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(StokesFlowSettings())
        self.settings.update(settings)

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
