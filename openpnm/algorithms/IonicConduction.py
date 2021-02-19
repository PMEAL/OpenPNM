import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.models.physics import generic_source_term as gst
from openpnm.utils import logging, Docorator, GenericSettings
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='IonicConductionSettings', sections=['Parameters'])
@docstr.dedent
class IonicConductionSettings(GenericSettings):
    r"""

    Parameters
    ----------
    phase : OpenPNM Phase object
        The phase on which the algorithm is to be run.
    quantity : str (default = ``'pore.mole_fraction'``)
          The name of the physical quantity to be calculated.
    conductance : str (default = ``'throat.ionic_conductance'``)
        The name of thepore-scale transport conductance values.  These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.
    charge_conservation : str (default = ``'electroneutrality'``)
        The assumption adopted to enforce charge conservation when
        performing ions transport simulations.

    Other Parameters
    ----------------

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """
    quantity = 'pore.potential'
    conductance = 'throat.ionic_conductance'
    charge_conservation = 'electroneutrality'
    cache_A = False
    cache_b = False


class IonicConduction(ReactiveTransport):
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

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings._update_settings_and_docs(IonicConductionSettings())
        self.settings.update(settings)

    @docstr.dedent
    def setup(self, phase=None, quantity='', conductance='',
              charge_conservation=None, **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Parameters
        ----------
        %(IonicConductionSettings.parameters)s

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
        Ps = (self['pore.all'] * np.isnan(self['pore.bc_value'])
              * np.isnan(self['pore.bc_rate']))
        mod = gst.charge_conservation
        phys = self.project.find_physics(phase=phase)
        phys[0].add_model(propname='pore.charge_conservation', model=mod,
                          phase=phase, p_alg=self, e_alg=e_alg,
                          assumption=self.settings['charge_conservation'])
        self.set_source(propname='pore.charge_conservation', pores=Ps)
