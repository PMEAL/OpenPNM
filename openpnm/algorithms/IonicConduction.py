import numpy as np
from openpnm.algorithms import ReactiveTransport
from openpnm.models.physics import generic_source_term as gst
from openpnm.utils import logging, Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.dedent
class IonicConductionSettings:
    r"""

    Parameters
    ----------
    %(ReactiveTransportSettings.parameters)s
    charge_conservation : str
        <??>
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
        self.settings = SettingsAttr(IonicConductionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)

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
