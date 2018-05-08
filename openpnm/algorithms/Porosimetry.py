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
                    'throat_entry_pressure': 'throat.capillary_pressure',
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

        throat_entry_pressure : string
            The dictionary key on the Phase object where the throat entry
            pressure values are stored.  The default is
            'throat.capillary_pressure'.  This is only accessed if the ``mode``
            is set to bond percolation.

        pore_volume : string
            The dictionary key containing the pore volume information. The
            default is 'pore.volume'.

        throat_volume : string
            The dictionary key containing the throat volume information.  The
            default is 'throat.volume'.

        late_pore_filling : string
            The name of the model used to determine partial pore filling as
            a function of applied pressure.

        late_throat_filling : string
            The name of the model used to determine partial throat filling as
            a function of applied pressure.

        """
        if phase:
            self.settings['phase'] = phase.name
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

    def run(self, points=25, start=None, stop=None):
        if self.settings['mode'] is not 'bond':
            raise Exception('Porosimetry must be run as bond percolation')
        if self.settings['access_limited'] is False:
            raise Exception('Porosimetry must be run as access limited')
        super().run(points=points, start=start, stop=stop)

    run.__doc__ = OrdinaryPercolation.run.__doc__

    def get_intrusion_data(self):
        r"""
        Obtain the numerical values of the calculated intrusion curve

        Returns
        -------
        A named-tuple containing arrays of applied capillary pressures and
        invading phase saturation.

        """
        net = self.project.network
        phase = self.project.find_phase(self)
        # Infer list of applied capillary pressures
        points = np.unique(self['throat.invasion_pressure'])
        # Add a low pressure point to the list to improve graph
        points = np.concatenate(([0], points))
        if points[-1] == np.inf:  # Remove infinity from PcPoints if present
            points = points[:-1]
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
            p_inv = self['pore.invasion_pressure'] <= p
            lpf = 1
            if self.settings['late_pore_filling']:
                # Set pressure on phase to current capillary pressure
                phase['pore.pressure'] = p
                # Regenerate corresponding physics model
                phase.regenerate_models(self.settings['late_pore_filling'])
                # Fetch partial filling fraction from phase object (0->1)
                lpf = phase[self.settings['late_pore_filling']]
            Vp = np.sum((Pvol*lpf)[p_inv])
            # Calculate filled throat volumes
            t_inv = self['throat.invasion_pressure'] <= p
            ltf = 1
            if self.settings['late_throat_filling']:
                # Set pressure on phase to current capillary pressure
                phase['throat.pressure'] = p
                # Regenerate corresponding physics model
                phase.regenerate_models(self.settings['late_pore_filling'])
                # Fetch partial filling fraction from phase object (0->1)
                ltf = phase[self.settings['late_pore_filling']]
            Vt = np.sum((Tvol*ltf)[t_inv])
            Vnwp_p.append(Vp)
            Vnwp_t.append(Vt)
            Vnwp_all.append(Vp + Vt)
        # Convert volumes to saturations by normalizing with total pore volume
        Snwp_all = [V/Total_vol for V in Vnwp_all]
        pc_curve = namedtuple('pc_curve', ('Pcap', 'Snwp'))
        data = pc_curve(points, Snwp_all)
        return data

    def plot_intrusion_curve(self):
        r"""
        Plot the percolation curve as the invader volume or number fraction vs
        the applied capillary pressure.

        """
        # Begin creating nicely formatted plot
        x, y = self.get_intrusion_data()
        fig = plt.figure()
        plt.semilogx(x, y, 'ko-')
        plt.ylabel('Invading Phase Saturation')
        plt.xlabel('Capillary Pressure')
        plt.grid(True)
        return fig
