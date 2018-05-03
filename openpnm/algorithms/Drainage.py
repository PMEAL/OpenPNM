import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from openpnm.algorithms import OrdinaryPercolation
from openpnm.core import logging
logger = logging.getLogger(__name__)


class Drainage(OrdinaryPercolation):
    r"""
    Simulates a capillary drainage experiment by applying a list of increasing
    capillary pressures and invading all throats that are accessible *and*
    invadable at the given pressure.

    Parameters
    ----------
    network : OpenPNM Network Object
        The network upon which the simulation will be run

    name : string (optional)
        The name to apply to the Algorithm for quick identification

    Notes
    -----
    This algorithm creates several arrays and adds them to its own dictionary.
    These include:

    **'pore(throat).inv_Pc'** : The applied capillary pressure at which each
    pore(throat) was invaded.

    **'pore(throat).inv_Seq'** : The sequence each pore(throat) was filled.

    Examples
    --------
    >>> import openpnm as op
    >>> pn = op.network.Cubic(shape=[20, 20, 20], spacing=10)
    >>> pn.add_boundary_pores(labels='top')
    >>> geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
    >>> water = op.phases.Water(network=pn)
    >>> air = op.phases.Air(network=pn)
    >>> phys = op.physics.GenericPhysics(network=pn, phase=water, geometry=geo)
    >>> phys.add_model(propname='throat.capillary_pressure',
    ...                model=op.models.physics.capillary_pressure.washburn)

    Once the basic Core objects are setup, the Algorithm can be created and
    and run as follows:

    >>> alg = op.algorithms.Drainage(network=pn)
    >>> alg.setup(phase=water)
    >>> alg.set_inlets(pores=pn.pores('top'))
    >>> alg.run()
    >>> data = alg.get_drainage_data()

    The ``data`` variable is a dictionary containing the numerical values of
    the resultant capillary pressure curve.  The available values can be
    inspected by typing ``data.keys()`` at the command line. These values can
    of course be plotted with matplotlib or exported to a graphing program of
    your choice.

    The final step is to utilize these results elsewhere in your OpenPNM
    simulation.  All Algorithms possess a method called ``return_results``
    which as the name suggests send the results to the correct locations.  For
    the 'Drainage Algorithm' this works as follows:

    >>> water.update(alg.results(Pc=5000))

    This command determines which pores and throats were filled at the applied
    capillary pressure of 5000, and creates 'pore.occupancy' and
    'throat.occupancy' arrays on the Phase object that was specfied as the
    'invading_phase' in the ``setup_method``.

    """

    def __init__(self, settings={}, **kwargs):
        super().__init__(**kwargs)
        # Apply class-specific default settings
        self.settings.update({'pore_volume': 'pore.volume',
                              'throat_volume': 'throat.volume',
                              'mode': 'bond',
                              'access_limited': True})
        # Override with any user specified settings
        self.settings.update(settings)

    def evaluate_late_pore_filling(self, Pc, Swp_init=0.75, eta=3.0,
                                   wetting_phase=False):
        r"""
        Compute the volume fraction of the phase in each pore given an initial
        wetting phase fraction (Swp_init) and a growth exponent (eta)
        returns the fraction of the pore volume occupied by wetting or
        non-wetting phase.
        Assumes Non-wetting phase displaces wetting phase
        """
        Swp = Swp_init*(self['pore.invasion_pressure']/Pc)**eta
        Swp[self['pore.invasion_pressure'] > Pc] = 1.0
        Snwp = 1-Swp
        if wetting_phase:
            return Swp
        else:
            return Snwp

    def _calc_fractional_filling(self, element, pressure):
        r"""
        Calculates the fractional filling of each pore or throat as the
        supplied capillary pressure

        Parameters
        ----------
        element : string
            Can either be 'pore' or 'throat' indicating which type of element
            to be calculated

        pressure : float
            The value of the capillary pressure to apply.  This value is sent
            the to Physics model stored in the ``pore_filling``
            attribute of the object.  This attribute is set during the call
            to ``get_drainage_data``.

        Notes
        -----
        The 'pore(throat)_filling' model must accept the applied capillary
        pressure as 'Pc'.  This is not customizable at the moment.
        """
        proj = self.project
        net = proj.network
        inv_phase = proj.phases()[self.settings['inv_phase']]
        if element == 'pore':
            key = self.settings['pore_filling']
            vol = self.settings['pore_volume']
        elif element == 'throat':
            key = self.settings['throat_filling']
            vol = self.settings['throat_volume']
        else:
            raise Exception('element must be either \'pore\' or \'throat\'')
        Snwp = np.zeros((self._count(element),))
        for phys in proj.find_physics(phase=inv_phase):
            # Update Physics model with current Pc
            temp_Pc = phys.models[key]['Pc']  # Store old Pc
            phys.models[key]['Pc'] = pressure
            # Regenerate Physics model and capture output locally
            Snwp[phys.Pnet] = phys.models[key].run()
            # Re-populate the residual element with the non-wetting phase
            if np.any(self[element+'.residual']):
                Snwp[self[element+'.residual']] = 1.0
            phys.models[key]['Pc'] = temp_Pc  # Return old Pc
        V = net[vol]*Snwp
        return V

    def get_drainage_data(self):
        r"""
        Obtain the numerical values of the calculate drainage pressure curve.

        Returns
        -------
        A named-tuple containing arrays of applied capillary pressures and
        non-wetting phase saturation.  A named-tuple means that the arrays
        can be accessed as named attributes like ``obj.Pcap``.

        """
        net = self.project.network
        # Infer list of applied capillary pressures
        PcPoints = np.unique(self['throat.invasion_pressure'])
        # Add a low pressure point to the list to improve graph
        PcPoints = np.concatenate(([0.9*PcPoints[0]], PcPoints))
        if PcPoints[-1] == np.inf:  # Remove infinity from PcPoints if present
            PcPoints = PcPoints[:-1]
        # Get pore and throat volumes
        Pvol = net[self.settings['pore_volume']]
        Tvol = net[self.settings['throat_volume']]
        Total_vol = np.sum(Pvol) + np.sum(Tvol)
        # Find cumulative filled volume at each applied capillary pressure
        Vnwp_t = []
        Vnwp_p = []
        Vnwp_all = []
        for Pc in PcPoints:
            # Calculate filled pore volumes
            p_inv = self['pore.invasion_pressure'] <= Pc
            if self.settings['pore_filling'] is None:
                Vp = np.sum(Pvol[p_inv])
            else:
                Vp = self._calc_fractional_filling(pressure=Pc,
                                                   element='pore')
                Vp = np.sum(Vp[p_inv])
            # Calculate filled throat volumes
            t_inv = self['throat.invasion_pressure'] <= Pc
            if self.settings['throat_filling'] is None:
                Vt = np.sum(Tvol[t_inv])
            else:
                Vt = self._calc_fractional_filling(pressure=Pc,
                                                   element='throat')
                Vt = np.sum(Vt[t_inv])
            Vnwp_p.append(Vp)
            Vnwp_t.append(Vt)
            Vnwp_all.append(Vp + Vt)
        # Convert volumes to saturations by normalizing with total pore volume
        Snwp_all = [V/Total_vol for V in Vnwp_all]
        pc_curve = namedtuple('pc_curve', ('Pcap', 'Snwp'))
        data = pc_curve(PcPoints, Snwp_all)
        return data

    def plot_drainage_curve(self):
        r"""
        Plot the drainage curve as the non-wetting phase saturation vs the
        applied capillary pressure.

        """
        # Begin creating nicely formatted plot
        data = self.get_drainage_data()
        xdata = data.Pcap
        ydata = data.Snwp
        fig = plt.figure()
        plt.semilogx(xdata, ydata, 'ko-')
        plt.ylabel('Invading Phase Saturation')
        plt.xlabel('Capillary Pressure')
        plt.grid(True)
        if np.amax(xdata) <= 1:
            plt.xlim(xmin=0, xmax=1)
        if np.amax(ydata) <= 1:
            plt.ylim(ymin=0, ymax=1)
        return fig
