import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from openpnm.algorithms import OrdinaryPercolation
from openpnm.core import logging
logger = logging.getLogger()

default_settings = {'pore_volume': 'pore.volume',
                    'throat_volume': 'throat.volume',
                    'mode': 'bond',
                    'access_limited': True,
                    'quantity': 'pressure',
                    'throat_entry_threshold': 'throat.entry_pressure',
                    'late_pore_filling': '',
                    'late_throat_filling': ''}


class Porosimetry(OrdinaryPercolation):

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        self.settings.update(default_settings)
        # Apply user settings, if any
        self.settings.update(settings)
        # Use the reset method to initialize all arrays
        self.reset()

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

        throat_entry_threshold : string
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
        if self.settings['mode'] is not 'bond':
            raise Exception('Porosimetry must be run as bond percolation')
        if self.settings['access_limited'] is False:
            raise Exception('Porosimetry must be run as access limited')
        super().run(points=points, start=start, stop=stop)

    run.__doc__ = OrdinaryPercolation.run.__doc__

    def get_intrusion_data(self, Pc=None):
        r"""
        Obtain the numerical values of the calculated intrusion curve

        Returns
        -------
        A named-tuple containing arrays of applied capillary pressures and
        invading phase saturation.

        """
        net = self.project.network
        if Pc is None:
            # Infer list of applied capillary pressures
            points = np.unique(self['throat.invasion_pressure'])
            # Add a low pressure point to the list to improve graph
            points = np.concatenate(([0], points))
            if points[-1] == np.inf:  # Remove infinity from points if present
                points = points[:-1]
        else:
            points = np.array(Pc)
        # Get pore and throat volumes
        Pvol = net[self.settings['pore_volume']]
        Tvol = net[self.settings['throat_volume']]
        Total_vol = np.sum(Pvol) + np.sum(Tvol)
        # Find cumulative filled volume at each applied capillary pressure
        Vnwp_t = []
        Vnwp_p = []
        Vnwp_all = []
        for p in points:
            # Calculate filled pore volumes
            p_inv, t_inv = self.results(p).values()
            Vp = np.sum(Pvol*p_inv)
            Vt = np.sum(Tvol*t_inv)
            Vnwp_p.append(Vp)
            Vnwp_t.append(Vt)
            Vnwp_all.append(Vp + Vt)
        # Convert volumes to saturations by normalizing with total pore volume
        Snwp_all = [V/Total_vol for V in Vnwp_all]
        pc_curve = namedtuple('pc_curve', ('Pcap', 'Snwp'))
        data = pc_curve(points, Snwp_all)
        return data

    def plot_intrusion_curve(self, fig=None):
        r"""
        Plot the percolation curve as the invader volume or number fraction vs
        the applied capillary pressure.

        """
        # Begin creating nicely formatted plot
        x, y = self.get_intrusion_data()
        if fig is None:
            fig = plt.figure()
        plt.semilogx(x, y, 'ko-')
        plt.ylabel('Invading Phase Saturation')
        plt.xlabel('Capillary Pressure')
        plt.grid(True)
        return fig

    def results(self, Pc):
        r"""
        """
        p_inv, t_inv = super().results(Pc).values()
        phase = self.project.find_phase(self)
        quantity = self.settings['quantity'].split('.')[-1]
        lpf = 1
        if self.settings['pore_partial_filling']:
            # Set pressure on phase to current capillary pressure
            phase['pore.'+quantity] = Pc
            # Regenerate corresponding physics model
            for phys in self.project.find_physics(phase=phase):
                phys.regenerate_models(self.settings['pore_partial_filling'])
            # Fetch partial filling fraction from phase object (0->1)
            lpf = phase[self.settings['pore_partial_filling']]
        # Calculate filled throat volumes
        ltf = 1
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
        return {'pore.occupancy': p_inv, 'throat.occupancy': t_inv}
