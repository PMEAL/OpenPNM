import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple
from openpnm.topotools import ispercolating
from openpnm.algorithms import OrdinaryPercolation


class SiteAndBondPercolation(OrdinaryPercolation):
    r"""
    Examples
    --------
    Start by importing all the necessary packages:

    >>> import openpnm as op
    >>> import scipy as sp
    >>> import matplotlib.pyplot as plt

    Initialize all the OpenPNM Objects and assign some basic properties:

    >>> shape = [50, 50]
    >>> pn = op.network.Cubic(shape=shape)
    >>> hg = op.phases.Mercury(network=pn)
    >>> phys = op.physics.GenericPhysics(network=pn, phase=hg, geometry=pn)
    >>> phys['pore.capillary_pressure'] = sp.rand(pn.Np)
    >>> phys['throat.capillary_pressure'] = sp.rand(pn.Nt)

    Initialize the algorithm, run setup, and generate a lit of invasion points:

    >>> mip = op.algorithms.SiteAndBondPercolation(network=pn)
    >>> mip.setup(phase=hg)
    >>> points = sp.linspace(0, 1, 20)

    Start by simulating access limited bond percolation:

    >>> mip.reset()
    >>> mip.settings['access_limited'] = True
    >>> mip.settings['mode'] = 'bond'
    >>> mip.set_inlets(pores=pn.pores(['left']))
    >>> mip.run(points)
    >>> data = mip.get_percolation_data()
    >>> fig = plt.subplot(2, 4, 1)
    >>> fig = plt.plot(*data, 'b-o')
    >>> fig = plt.subplot(2, 4, 2)
    >>> fig = plt.imshow(sp.reshape(mip['pore.invasion_pressure'], shape))
    >>> fig = plt.title('Access Limited Bond Percolation')

    Clear the previous results with reset and rerun without access limitations:

    >>> mip.reset()
    >>> mip.settings['access_limited'] = False
    >>> mip.settings['mode'] = 'bond'
    >>> mip.run(points)
    >>> data = mip.get_percolation_data()
    >>> fig = plt.subplot(2, 4, 3)
    >>> fig = plt.plot(*data, 'r-o')
    >>> fig = plt.subplot(2, 4, 4)
    >>> fig = plt.imshow(sp.reshape(mip['pore.invasion_pressure'], shape))
    >>> fig = plt.title('Normal Bond Percolation')

    Now try site percolation with access limitations:

    >>> mip.reset()
    >>> mip.settings['access_limited'] = True
    >>> mip.settings['mode'] = 'site'
    >>> mip.set_inlets(pores=pn.pores(['left']))
    >>> mip.run(points)
    >>> data = mip.get_percolation_data()
    >>> fig = plt.subplot(2, 4, 5)
    >>> fig = plt.plot(*data, 'g-o')
    >>> fig = plt.subplot(2, 4, 6)
    >>> fig = plt.imshow(sp.reshape(mip['pore.invasion_pressure'], shape))
    >>> fig = plt.title('Access Limited Site Percolation')

    And finally, site percolation without access limitations:

    >>> mip.reset()
    >>> mip.settings['access_limited'] = False
    >>> mip.settings['mode'] = 'site'
    >>> mip.run(points)
    >>> data = mip.get_percolation_data()
    >>> fig = plt.subplot(2, 4, 7)
    >>> fig = plt.plot(*data, 'c-o')
    >>> fig = plt.subplot(2, 4, 8)
    >>> fig = plt.imshow(sp.reshape(mip['pore.invasion_pressure'], shape))
    >>> fig = plt.title('Normal Site Percolation')
    """

    def get_percolation_threshold(self):
        r"""
        """
        if np.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores must be specified first')
        if np.sum(self['pore.outlets']) == 0:
            raise Exception('Outlet pores must be specified first')
        else:
            Pout = self['pore.outlets']
        # Do a simple check of pressures on the outlet pores first...
        if self.settings['access_limited']:
            thresh = np.amin(self['pore.invasion_pressure'][Pout])
        else:
            raise Exception('This is currently only implemented for access ' +
                            'limited simulations')
        return thresh

    def is_percolating(self, applied_pressure):
        r"""
        Returns a True or False value to indicate if a percolating cluster
        spans between the inlet and outlet pores that were specified.

        Parameters
        ----------
        applied_pressure : scalar, float
            The pressure at which percolation should be checked

        Returns
        -------
        A simple boolean True or False if percolation has occured or not.

        """
        if np.sum(self['pore.inlets']) == 0:
            raise Exception('Inlet pores must be specified first')
        else:
            Pin = self['pore.inlets']
        if np.sum(self['pore.outlets']) == 0:
            raise Exception('Outlet pores must be specified first')
        else:
            Pout = self['pore.outlets']
        # Do a simple check of pressures on the outlet pores first...
        if np.amin(self['pore.invasion_pressure'][Pout]) > applied_pressure:
            val = False
        else:  # ... and do a rigorous check only if necessary
            mask = self['throat.invasion_pressure'] < applied_pressure
            am = self.project.network.create_adjacency_matrix(weights=mask,
                                                              fmt='coo')
            val = ispercolating(am=am, mode=self.settings['mode'],
                                inlets=Pin, outlets=Pout)
        return val

    def get_percolation_data(self):
        r"""
        Obtain the numerical values of the calculated percolation curve

        Returns
        -------
        A named-tuple containing arrays of applied capillary pressures and
        invading phase saturation.

        """
        net = self.project.network
        # Infer list of applied capillary pressures
        points = np.unique(self['throat.invasion_pressure'])
        # Add a low pressure point to the list to improve graph
        points = np.concatenate(([0], points))
        if points[-1] == np.inf:  # Remove infinity from PcPoints if present
            points = points[:-1]
        # Get pore and throat volumes
        if self.settings['pore_volume']:
            Pvol = net[self.settings['pore_volume']]
        else:
            Pvol = np.ones(shape=(self.Np, ), dtype=int)
        if self.settings['throat_volume']:
            Tvol = net[self.settings['throat_volume']]
        else:
            Tvol = np.zeros(shape=(self.Nt, ), dtype=int)
        Total_vol = np.sum(Pvol) + np.sum(Tvol)
        # Find cumulative filled volume at each applied capillary pressure
        Vnwp_t = []
        Vnwp_p = []
        Vnwp_all = []
        for p in points:
            # Calculate filled pore volumes
            p_inv = self['pore.invasion_pressure'] <= p
            Vp = np.sum(Pvol[p_inv])
            # Calculate filled throat volumes
            t_inv = self['throat.invasion_pressure'] <= p
            Vt = np.sum(Tvol[t_inv])
            Vnwp_p.append(Vp)
            Vnwp_t.append(Vt)
            Vnwp_all.append(Vp + Vt)
        # Convert volumes to saturations by normalizing with total pore volume
        Snwp_all = [V/Total_vol for V in Vnwp_all]
        pc_curve = namedtuple('pc_curve', ('Pcap', 'Snwp'))
        data = pc_curve(points, Snwp_all)
        return data

    def plot_percolation_curve(self):
        r"""
        Plot the percolation curve as the invader volume or number fraction vs
        the applied capillary pressure.

        """
        # Begin creating nicely formatted plot
        data = self.get_percolation_data()
        xdata = data.Pcap
        ydata = data.Snwp
        fig = plt.figure()
        plt.plot(xdata, ydata, 'ko-')
        plt.ylabel('Invading Phase Saturation')
        plt.xlabel('Capillary Pressure')
        plt.grid(True)
        if np.amax(xdata) <= 1:
            plt.xlim(xmin=0, xmax=1)
        if np.amax(ydata) <= 1:
            plt.ylim(ymin=0, ymax=1)
        return fig
