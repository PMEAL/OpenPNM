# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationBasic2: Simple IP
===============================================================================


"""

import scipy as sp
import heapq as hq
from OpenPNM.Algorithms import GenericAlgorithm
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class InvasionPercolation2(GenericAlgorithm):
    r"""
    A classic/basic invasion percolation algorithm.  This version of the
    algorithm allows for incremental filling of a specified number of throats,
    and can be called repeatedly.  This is useful for inside a for loop where
    some property is calculated at each invasion level.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network upon which the invasion should occur.

    Usage
    -----
    After instantiating an IP object, the next step is to call the ``setup``
    method to specify the various parameters such as the inlet pores.  To
    finally execute the algorithm, use the ``run`` method.  See the docstrings
    of each method for more details.

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._tcount = 0
        self['throat.invaded'] = -sp.ones_like(self._net.Ts)
        self['pore.invaded'] = -sp.ones_like(self._net.Ps)

    def setup(self, phase, p_inlets, throat_prop='throat.capillary_pressure'):
        r"""
        Specify the overall parameters for the algorithm

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        p_inlets : array_like
            The list of inlet pores from which the Phase can enter the Network.

        throat_prop : string
            The name of the throat property containing the capillary entry
            pressure.  The default is 'throat.capillary_pressure'.

        Notes
        -----
        When ``setup`` is called the entire Algorithm is reset, so this can be
        used to repeat a calculation if desired.
        """
        p_inlets = sp.array(p_inlets, ndmin=1)
        self._phase = phase
        # Setup arrays and info
        t_entry = phase[throat_prop]
        # Indices into t_entry giving a sorted list
        self['throat.sorted'] = sp.argsort(t_entry, axis=0)
        self['throat.order'] = sp.zeros_like(self['throat.sorted'])
        self['throat.order'][self['throat.sorted']] = sp.arange(0, self._net.Nt)

        self['pore.invaded'][p_inlets] = self._tcount
        # Perform initial analysis on input pores
        Ts = self._net.find_neighbor_throats(pores=p_inlets, mode='intersection')
        # Set throats connecting inlet pores to filled
        self['throat.invaded'][Ts] = self._tcount
        # Set other throats as potential invaded thorats
        Ts = self._net.find_neighbor_throats(pores=p_inlets, mode='not_intersection')
        # Add throats to the queue
        self._queue = []
        [hq.heappush(self._queue, T) for T in self['throat.order'][Ts]]
        hq.heapify(self._queue)

    def run(self, nsteps=None):
        r"""
        Perform the algorithm

        Parameters
        ----------
        nsteps : int
            The number of throats to invaded.  By default, all throats will be
            invaded.

        Returns
        -------
        filled : boolean
            This method returns False if there are throats remaining to be
            invaded, and True once all throats have been invaded.  This can be
            used for controlling a while loop to successively invade the
            Network.

        Examples
        --------
        To invade the Network in successive steps, simply call the run method
        inside a ``while`` loop, using the returned ``filled`` value to halt
        when complete:

        .. code-block:: python

            filled = False
            while not filled:
                filled = IP.run(nsteps=10)
                # Perform some other calculations here...

        """
        if nsteps is None:
            nsteps = sp.inf

        # Begin looping through accessible throats
        old_count = self._tcount
        while len(self._queue) > 0 and (self._tcount < (old_count + nsteps)):
            self._tcount += 1
            # Find throat at the top of the queue
            t = hq.heappop(self._queue)
            # Extract actual throat number
            t_next = self['throat.sorted'][t]
            # If throat is duplicated
            while len(self._queue) > 0 and self._queue[0] == t:
                # Note: Preventing duplicate entries below might save some time here
                t = hq.heappop(self._queue)
            self['throat.invaded'][t_next] = self._tcount
            # Find pores connected to newly invaded throat
            Ps = self._net['throat.conns'][t_next]
            # Remove already invaded pores from Ps
            Ps = Ps[self['pore.invaded'][Ps] < 0]
            if len(Ps) > 0:
                self['pore.invaded'][Ps] = self._tcount
                Ts = self._net.find_neighbor_throats(pores=Ps)  # Find connected throats
                Ts = Ts[self['throat.invaded'][Ts] < 0]  # Remove already invaded throats from Ts
                [hq.heappush(self._queue, T) for T in self['throat.order'][Ts]]
        complete = len(self._queue) == 0
        return complete

    def return_results(self, pores=[], throats=[]):
        r"""
        Places the results of the IP simulation into the Phase object.

        Parameters
        ----------
        pores and throats : array_like
            The list of pores and throats whose values should be returned to
            the Phase object.  Default is all of them.

        Returns
        -------
        invasion_sequence : array_like
            The sequence in which each pore and throat is invaded  This depends
            on the inlet locations.  All inlets are invaded at step 0.  It is
            possible to recontruct an animation of the invasion process, in
            Paraview for instance, using this sequence information.

        """
        pores = sp.array(pores, ndmin=1)
        throats = sp.array(throats, ndmin=1)
        if len(pores) == 0:
            pores = self.Ps
        if len(throats) == 0:
            throats = self.Ts
        self._phase['throat.invasion_sequence'] = sp.nan
        self._phase['pore.invasion_sequence'] = sp.nan
        self._phase['throat.invasion_sequence'][throats] = \
            self['throat.invaded'][throats]
        self._phase['pore.invasion_sequence'][pores] = \
            self['pore.invaded'][pores]

    def apply_flow(self, flowrate):
        r"""
        Convert the invaded sequence into an invaded time for a given flow rate
        considering the volume of invaded pores and throats.

        Parameters
        ----------
        flowrate : float
            The flow rate of the injected fluid

        Returns
        -------
        Creates a throat array called 'invasion_time' in the Algorithm
        dictionary

        """
        P12 = self._net['throat.conns']
        a = self['throat.invaded']
        b = sp.argsort(self['throat.invaded'])
        P12_inv = self['pore.invaded'][P12]
        # Find if the connected pores were invaded with or before each throat
        P1_inv = P12_inv[:, 0] == a
        P2_inv = P12_inv[:, 1] == a
        c = sp.column_stack((P1_inv, P2_inv))
        d = sp.sum(c, axis=1, dtype=bool)  # List of Pores invaded with each throat
        # Find volume of these pores
        P12_vol = sp.zeros((self.Nt,))
        P12_vol[d] = self._net['pore.volume'][P12[c]]
        # Add invaded throat volume to pore volume (if invaded)
        T_vol = P12_vol + self._net['throat.volume']
        # Cumulative sum on the sorted throats gives cumulated inject volume
        e = sp.cumsum(T_vol[b] / flowrate)
        t = sp.zeros((self.Nt,))
        t[b] = e  # Convert back to original order
        self._phase['throat.invasion_time'] = t
