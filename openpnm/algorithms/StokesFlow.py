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
        self.settings.update({'quantity': 'pore.pressure',
                              'conductance': 'throat.hydraulic_conductance'})
        self.settings.update(settings)
        super().__init__(**kwargs)

    def setup(self, phase=None, quantity='', conductance='', **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase on which the algorithm is to be run.  If no value is
            given, the existing value is kept.

        quantity : string
            The name of the physical quantity to be calcualted.  If no value is
            given, the existing value is kept.  The default value is
            ``'pore.pressure'``.

        conductance : string
            The name of the pore-scale transport conductance values.  These
            are typically calculate by a model attached to a *Physics* object
            associated with the given *Phase*.  If no value is given, the
            existing value is kept.  The default value is
            ``'throat.hydraulic_conductance'``.

        Notes
        -----
        Any additional arguments are added to the ``settings`` dictionary of
        the object.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity': quantity]
        if conductance:
            self.settings['conductance': conductance]
        super().setup(**kwargs)

    def calc_eff_permeability(self):
        r"""
        This calculates the effective permeability in this linear transport
        algorithm.
        """
        phase = self.project.phases()[self.settings['phase']]
        d_normal = self._calc_eff_prop()
        self._eff_property = d_normal * sp.mean(phase['pore.viscosity'])
        return self._eff_property
