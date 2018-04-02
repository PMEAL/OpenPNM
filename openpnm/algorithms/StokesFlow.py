import scipy as sp
from openpnm.algorithms import ReactiveTransport
from openpnm.core import logging
logger = logging.getLogger(__name__)


class StokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.  The 2
    main roles of this subclass are to set the default property names and to
    implement a method for calculating the hydraulic permeability of the
    network.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update({'quantity': 'pore.pressure',
                              'conductance': 'throat.hydraulic_conductance'})
        self.settings.update(settings)

    def calc_eff_permeability(self):
        r"""
        This calculates the effective permeability in this linear
        transport algorithm.
        """
        phase = self.project.phases[self.settings['phase']]
        d_normal = self._calc_eff_prop()
        self._eff_property = d_normal / sp.mean(phase['pore.viscosity'])
        return self._eff_property
