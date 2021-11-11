import numpy as np
from openpnm.algorithms import OrdinaryPercolation
from openpnm.utils import logging, SettingsAttr, Docorator
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.dedent
class PorosimetrySettings:
    r"""
    %(OrdinaryPercolationSettings.parameters)s
    pore_partial_filling : string
        The name of the model used to determine partial pore filling as
        a function of applied pressure.
    throat_partial_filling : string
        The name of the model used to determine partial throat filling as
        a function of applied pressure.
    """
    quantity = 'pore.pressure'
    pore_partial_filling = ''
    throat_partial_filling = ''


class Porosimetry(OrdinaryPercolation):
    r"""
    Simulates mercury instrustion porosimetry using ordinary percolation

    Parameters
    ----------
    network : OpenPNM Network object
        The Network upon which this simulation should be run
    name : string, optional
        An identifying name for the object.  If none is given then one is
        generated.

    Notes
    -----
    Mercury intrusion progresses by applying increasing pressures to the
    invading mercury phase, and measuring the resultant volume of invading
    fluid.  This corresponds directly to an ordinary percolation process,
    with access limitations enabled.

    See Also
    --------
    OrdinaryPercolation


    """

    def __init__(self, settings={}, phase=None, **kwargs):
        self.settings = SettingsAttr(PorosimetrySettings, settings)
        super().__init__(settings=self.settings, **kwargs)
        # Use the reset method to initialize all arrays
        self.reset()
        if phase is not None:
            self.settings['phase'] = phase.name

    def set_partial_filling(self, propname):
        r"""
        Define which pore filling model to apply.

        Parameters
        ----------
        propname : string
            Dictionary key on the physics object(s) containing the pore
            filling model(s) to apply.

        Notes
        -----
        It is assumed that these models are functions of the `quantity`
        specified in the algorithms settings.  This values is applied to the
        corresponding phase just prior to regenerating the given pore-scale
        model(s).

        """
        if propname.startswith('pore'):
            self.settings['pore_partial_filling'] = propname
        if propname.startswith('throat'):
            self.settings['throat_partial_filling'] = propname

    def run(self, points=25, start=None, stop=None):
        if self.settings['mode'] != 'bond':
            raise Exception('Porosimetry must be run as bond percolation')
        if self.settings['access_limited'] is False:
            raise Exception('Porosimetry must be run as access limited')
        super().run(points=points, start=start, stop=stop)

    run.__doc__ = OrdinaryPercolation.run.__doc__

    def results(self, Pc=None):
        r"""
        """
        if Pc is None:
            p_inv = self['pore.invasion_pressure']
            t_inv = self['throat.invasion_pressure']
            results = {'pore.invasion_pressure': p_inv,
                       'throat.invasion_pressure': t_inv}
        else:
            p_inv, t_inv = super().results(Pc).values()
            phase = self.project[self.settings.phase]
            quantity = self.settings['quantity'].split('.')[-1]
            lpf = np.array([1])
            if self.settings['pore_partial_filling']:
                # Set pressure on phase to current capillary pressure
                phase['pore.'+quantity] = Pc
                # Regenerate corresponding physics model
                for phys in self.project.find_physics(phase=phase):
                    phys.regenerate_models(self.settings['pore_partial_filling'])
                # Fetch partial filling fraction from phase object (0->1)
                lpf = phase[self.settings['pore_partial_filling']]
            # Calculate filled throat volumes
            ltf = np.array([1])
            if self.settings['throat_partial_filling']:
                # Set pressure on phase to current capillary pressure
                phase['throat.'+quantity] = Pc
                # Regenerate corresponding physics model
                for phys in self.project.find_physics(phase=phase):
                    phys.regenerate_models(self.settings['throat_partial_filling'])
                # Fetch partial filling fraction from phase object (0->1)
                ltf = phase[self.settings['throat_partial_filling']]
            p_inv = p_inv*lpf
            t_inv = t_inv*ltf
            results = {'pore.occupancy': p_inv, 'throat.occupancy': t_inv}
        return results
