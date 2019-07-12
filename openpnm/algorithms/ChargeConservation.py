import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.models.physics import generic_source_term as gst
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class ChargeConservation(ReactiveTransport):
    r"""
    A class to enforce charge conservation in ionic transport simulations.

    Parameters
    ----------
    network : OpenPNM Network object
        The network on which this algorithm operates

    project : OpenPNM Project object
        Either a network or a project must be specified

    name : string, optional
        A unique name to give the object for easier identification.  If not
        given, one is generated.
    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.potential',
                   'conductance': 'throat.ionic_conductance',
                   'charge_conservation': 'electroneutrality',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': '',
                                            'charge_conservation': ''},
                           'set_rate_BC':  {'pores': None,
                                            'values': None},
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase=None, quantity='', conductance='',
              charge_conservation=None, **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase on which the algorithm is to be run.

        quantity : string
            (default is ``'pore.mole_fraction'``)  The name of the physical
            quantity to be calculated.

        conductance : string
            (default is ``'throat.diffusive_conductance'``) The name of the
            pore-scale transport conductance values.  These are typically
            calculated by a model attached to a *Physics* object associated
            with the given *Phase*.

        charge_conservation : string
            The assumption adopted to enforce charge conservation when
            performing ions transport simulations (default is
            "electroneutrality").

        Notes
        -----
        Any additional arguments are added to the ``settings`` dictionary of
        the object.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        if charge_conservation:
            self.settings['charge_conservation'] = charge_conservation
        super().setup(**kwargs)

    def _charge_conservation_eq_source_term(self, e_alg):
        # Source term for Poisson or charge conservation (electroneutrality) eq
        phase = self.project.phases()[self.settings['phase']]
        Ps = (self['pore.all'] * np.isnan(self['pore.bc_value']) *
              np.isnan(self['pore.bc_rate']))
        mod = gst.charge_conservation
        phys = self.project.find_physics(phase=phase)
        phys[0].add_model(propname='pore.charge_conservation', model=mod,
                          phase=phase, p_alg=self, e_alg=e_alg,
                          assumption=self.settings['charge_conservation'])
        self.set_source(propname='pore.charge_conservation', pores=Ps)
