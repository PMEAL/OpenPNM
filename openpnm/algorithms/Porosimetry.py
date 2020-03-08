from openpnm.algorithms import OrdinaryPercolation
from openpnm.utils import logging
import numpy as np
logger = logging.getLogger(__name__)


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

    project : OpenPNM Project object
        Either a Network or a Project must be specified

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
        def_set = {'phase': None,
                   'pore_volume': 'pore.volume',
                   'throat_volume': 'throat.volume',
                   'mode': 'bond',
                   'access_limited': True,
                   'quantity': 'pressure',
                   'throat_entry_pressure': 'throat.entry_pressure',
                   'pore_volume': 'pore.volume',
                   'throat_volume': 'throat.volume',
                   'late_pore_filling': '',
                   'late_throat_filling': '',
                   'gui': {'setup': {'phase': None,
                                     'quantity': '',
                                     'throat_entry_pressure': '',
                                     'pore_volume': '',
                                     'throat_volume': '',
                                     'late_pore_filling': '',
                                     'late_throat_filling': ''},
                           'set_inlets':   {'pores': None,
                                            'overwrite': False},
                           'set_outlets':  {'pores': None,
                                            'overwrite': False},
                           'set_residual': {'pores': None,
                                            'throats': None,
                                            'overwrite': False}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        # Apply user settings, if any
        self.settings.update(settings)
        # Use the reset method to initialize all arrays
        self.reset()
        if phase is not None:
            self.setup(phase=phase)

    def setup(self,
              phase=None,
              quantity='',
              throat_entry_pressure='',
              pore_volume='',
              throat_volume='',
              late_pore_filling='',
              late_throat_filling=''):
        r"""
        Used to specify necessary arguments to the simulation.  This method is
        useful for resetting the algorithm or applying more explicit control.

        Parameters
        ----------
        phase : OpenPNM Phase object
            The Phase object containing the physical properties of the invading
            fluid.

        quantity : string
            The name of the quantity calculated by this algorithm.  This is
            used for instance, by the late pore and throat filling models
            to indicate the prevailing fluid pressure in the invading phase
            for calculating the extent of filling.  The default is
            'pressure'.  Note that there is no need to specify 'pore' and/or
            'throat' with this as the given value will apply to both.

        throat_entry_pressure : string
            The dictionary key on the Phase object where the throat entry
            pressure values are stored.  The default is
            'throat.entry_pressure'.

        pore_volume : string
            The dictionary key containing the pore volume information. The
            default is 'pore.volume'.

        throat_volume : string
            The dictionary key containing the throat volume information.  The
            default is 'throat.volume'.

        pore_partial_filling : string
            The name of the model used to determine partial pore filling as
            a function of applied pressure.

        throat_partial_filling : string
            The name of the model used to determine partial throat filling as
            a function of applied pressure.

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if throat_entry_pressure:
            self.settings['throat_entry_pressure'] = throat_entry_pressure
            phase = self.project.find_phase(self)
            self['throat.entry_pressure'] = phase[throat_entry_pressure]
        if pore_volume:
            self.settings['pore_volume'] = pore_volume
        if throat_volume:
            self.settings['throat_volume'] = throat_volume
        if late_pore_filling:
            self.settings['late_pore_filling'] = late_pore_filling
        if late_throat_filling:
            self.settings['late_throat_filling'] = late_throat_filling

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
            phase = self.project.find_phase(self)
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
