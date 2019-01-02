import heapq as hq
import scipy as sp
import numpy as np
from openpnm.algorithms import GenericAlgorithm
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class InvasionPercolation(GenericAlgorithm):
    r"""
    A classic/basic invasion percolation algorithm optimized for speed.

    Parameters
    ----------
    network : OpenPNM Network object
        The Network upon which the invasion will occur.

    Notes
    ----
    This algorithm uses a binary heap to store all a list of all accessible
    throats, sorted according to entry pressure.  This means that item [0] in
    the heap is the most easily invaded throat, so looking up which throat
    to invade next is computationally trivial.  In order to keep the list
    sorted new throats to the list takes more time, however, the heap data
    structure is very efficient at this.  Interested users can consult the
    wikipedia page on `binary heaps
    <https://en.wikipedia.org/wiki/Binary_heap>`_ for more information.


    Examples
    --------
    Start by importing the usual packages:

    >>> import openpnm as op
    >>> import scipy as sp
    >>> import matplotlib.pyplot as plt

    Create 2D cubic network for easier visualizaiton:

    >>> S = sp.array([100, 100, 1])
    >>> pn = op.network.Cubic(shape=S, spacing=0.0001, name='pn11')

    Add a basic geometry:

    >>> geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

    Create an invading phase, and attach the capillary pressure model:

    >>> water = op.phases.Water(network=pn)
    >>> water.add_model(propname='throat.entry_pressure',
    ...                 model=op.models.physics.capillary_pressure.washburn)

    Initialize an invasion percolation object and define inlets:

    >>> ip = op.algorithms.InvasionPercolation(network=pn)
    >>> ip.setup(phase=water)
    >>> ip.set_inlets(pores=0)
    >>> ip.run()

    After running the algorithm the invading phase configuration at a given
    saturation can be obtained and assigned to the phase object:

    >>> water.update(ip.results(Snwp=0.5))

    Because it was a 2D network it's easy to quickly visualize the invasion
    pattern as an image for verification:

    .. note::

        Because the network is 2D and cubic, an image can be generated with
        color corresponding to a value.  The following plots the entire
        invasion sequence, and the water configuraiton at Snwp = 0.5.

        ``plt.subplot(1, 2, 1)``

        ``plt.imshow(sp.reshape(ip['pore.invasion_sequence'], newshape=S[S > 1]))``

        ``plt.subplot(1, 2, 2)``

        ``plt.imshow(sp.reshape(water['pore.occupancy'], newshape=S[S > 1]))``

    """
    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'pore_volume': 'pore.volume',
                   'throat_volume': 'throat.volume',
                   'entry_pressure': 'throat.entry_pressure',
                   'gui': {'setup':          {'phase': None,
                                              'entry_pressure': '',
                                              'pore_volume': '',
                                              'throat_volume': ''},
                           'set_inlets':     {'pores': None,
                                              'overwrite': False},
                           'apply_trapping': {'outlets': None}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    def setup(self, phase, entry_pressure='', pore_volume='', throat_volume=''):
        r"""
        Set up the required parameters for the algorithm

        Parameters
        ----------
        phase : OpenPNM Phase object
            The phase to be injected into the Network.  The Phase must have the
            capillary entry pressure values for the system.

        entry_pressure : string
            The dictionary key to the capillary entry pressure.  If none is
            supplied then the current value is retained. The default is
            'throat.capillary_pressure'.

        pore_volume : string
            The dictionary key to the pore volume.  If none is supplied then
            the current value is retained. The default is 'pore.volume'.

        throat_volume : string
            The dictionary key to the throat volume.  If none is supplied then
            the current value is retained. The default is 'throat.volume'.

        """
        self.settings['phase'] = phase.name
        if pore_volume:
            self.settings['pore_volume'] = pore_volume
        if throat_volume:
            self.settings['throat_volume'] = throat_volume
        if entry_pressure:
            self.settings['entry_pressure'] = entry_pressure

        # Setup arrays and info
        self['throat.entry_pressure'] = phase[self.settings['entry_pressure']]
        # Indices into t_entry giving a sorted list
        self['throat.sorted'] = sp.argsort(self['throat.entry_pressure'], axis=0)
        self['throat.order'] = 0
        self['throat.order'][self['throat.sorted']] = sp.arange(0, self.Nt)
        self['throat.invasion_sequence'] = -1
        self['pore.invasion_sequence'] = -1
        self._tcount = 0

    def set_inlets(self, pores=[], overwrite=False):
        r"""

        Parameters
        ----------
        pores : array_like
            The list of inlet pores from which the Phase can enter the Network
        """
        if overwrite:
            self['pore.invasion_sequence'] = -1
        self['pore.invasion_sequence'][pores] = 0

        # Perform initial analysis on input pores
        Ts = self.project.network.find_neighbor_throats(pores=pores)
        self.queue = []
        [hq.heappush(self.queue, T) for T in self['throat.order'][Ts]]

    def run(self, n_steps=None):
        r"""
        Perform the algorithm

        Parameters
        ----------
        n_steps : int
            The number of throats to invaded during this step

        """
        if n_steps is None:
            n_steps = sp.inf

        queue = self.queue
        if len(queue) == 0:
            logger.warn('queue is empty, this network is fully invaded')
            return
        t_sorted = self['throat.sorted']
        t_order = self['throat.order']
        t_inv = self['throat.invasion_sequence']
        p_inv = self['pore.invasion_sequence']

        count = 0
        while (len(queue) > 0) and (count < n_steps):
            # Find throat at the top of the queue
            t = hq.heappop(queue)
            # Extract actual throat number
            t_next = t_sorted[t]
            t_inv[t_next] = self._tcount
            # If throat is duplicated
            while len(queue) > 0 and queue[0] == t:
                # Note: Preventing duplicate entries below might save some time
                t = hq.heappop(queue)
            # Find pores connected to newly invaded throat
            Ps = self.project.network['throat.conns'][t_next]
            # Remove already invaded pores from Ps
            Ps = Ps[p_inv[Ps] < 0]
            if len(Ps) > 0:
                p_inv[Ps] = self._tcount
                Ts = self.project.network.find_neighbor_throats(pores=Ps)
                Ts = Ts[t_inv[Ts] < 0]  # Remove invaded throats from Ts
                [hq.heappush(queue, T) for T in t_order[Ts]]
            count += 1
            self._tcount += 1
        self['throat.invasion_sequence'] = t_inv
        self['pore.invasion_sequence'] = p_inv

    def results(self, Snwp):
        r"""
        Returns the phase configuration at the specified non-wetting phase
        (invading phase) saturation.

        Parameters
        ----------
        Snwp : scalar, between 0 and 1
            The network saturation for which the phase configuration is
            desired.

        Returns
        -------
        Two dictionary containing arrays that describe the pore and throat
        distribution at the given saturation.  Specifically, these are:

        **'pore.occupancy'** : 1 indicates the pores is invaded and 0
        otherwise.

        **'throat.occupancy'** : Same as described above but for throats.

        """
        net = self.project.network
        P12 = net['throat.conns']
        # Fetch void volume for pores and throats
        Vp = net[self.settings['pore_volume']]
        Vt = net[self.settings['throat_volume']]
        # Fetch the order of filling
        Np = self['pore.invasion_sequence']
        Nt = self['throat.invasion_sequence']
        # Create Nt-long mask of which pores were filled when throat was filled
        Pinv = (Np[P12].T == Nt).T
        # If a pore and throat filled together, find combined volume
        Vinv = sp.vstack(((Pinv*Vp[P12]).T, Vt)).T
        Vinv = sp.sum(Vinv, axis=1)
        # Convert to cumulative volume filled as each throat is invaded
        x = sp.argsort(Nt)  # Find order throats were invaded
        Vinv_cum = np.cumsum(Vinv[x])
        # Normalized cumulative volume filled into saturation
        S = Vinv_cum/(Vp.sum() + Vt.sum())
        # Find throat invasion step where Snwp was reached
        try:
            N = sp.where(S < Snwp)[0][-1]
        except:
            N = -np.inf
        data = {'pore.occupancy': Np <= N, 'throat.occupancy': Nt <= N}
        return data

    def apply_trapping(self, outlets):
        """
        Apply trapping based on algorithm described by Y. Masson [1].
        It is applied as a post-process and runs the percolation algorithm in
        reverse assessing the occupancy of pore neighbors. Consider the
        following scenario when running standard IP without trapping,
        3 situations can happen after each invasion step:

            * The number of defending clusters stays the same and clusters can
            shrink
            * A cluster of size one is suppressed
            * A cluster is split into multiple clusters

        In reverse the following opposite situations can happen:

            * The number of defending clusters stays the same and clusters can
            grow
            * A cluster of size one is created
            * Mutliple clusters merge into one cluster

        With trapping the reversed rules are adjusted so that only clusters
        that do not connect to a sink can grow and merge. At the point that a
        neighbor connected to a sink is touched the trapped cluster stops
        growing as this is the point of trapping in forward invasion time.

        Logger info displays the invasion sequence and pore index and a message
        with condition number based on the modified trapping rules and the
        assignment of the pore to a given cluster.

        Initially all invaded pores are given cluster label -1
        Outlets / Sinks are given -2
        New clusters that grow into fully trapped clusters are either
        identified at the point of breakthrough or grow from nothing if the
        full invasion sequence is run, they are assigned numbers from 0 up.

        Ref:
        [1] Masson, Y., 2016. A fast two-step algorithm for invasion
        percolation with trapping. Computers & Geosciences, 90, pp.41-48

        Parameters
        ----------
        outlets : list or array of pore indices for defending fluid to escape
        through

        Returns
        -------
        Creates a throat array called 'pore.clusters' in the Algorithm
        dictionary. Any positive number is a trapped cluster

        Also creates 2 boolean arrays Np and Nt long called '<element>.trapped'
        """
        # First see if network is fully invaded
        net = self.project.network
        invaded_ps = self['pore.invasion_sequence'] > -1
        if ~np.all(invaded_ps):
            # Put defending phase into clusters
            clusters = net.find_clusters2(~invaded_ps)
            # Identify clusters that are connected to an outlet and set to -2
            # -1 is the invaded fluid
            # -2 is the defender fluid able to escape
            # All others now trapped clusters which grow as invasion is reversed
            out_clusters = sp.unique(clusters[outlets])
            for c in out_clusters:
                if c >= 0:
                    clusters[clusters == c] = -2
        else:
            # Go from end
            clusters = np.ones(net.Np, dtype=int)*-1
            clusters[outlets] = -2

        # Turn into a list for indexing
        inv_seq = np.vstack((self['pore.invasion_sequence'].astype(int),
                             np.arange(0, net.Np, dtype=int))).T
        # Reverse sort list
        inv_seq = inv_seq[inv_seq[:, 0].argsort()][::-1]
        next_cluster_num = np.max(clusters)+1
        # For all the steps after the inlets are set up to break-through
        # Reverse the sequence and assess the neighbors cluster state
        stopped_clusters = np.zeros(net.Np, dtype=bool)
        all_neighbors = net.find_neighbor_pores(net.pores(), flatten=False,
                                                include_input=True)
        for un_seq, pore in inv_seq:
            if pore not in outlets and un_seq > 0:  # Skip inlets and outlets
                nc = clusters[all_neighbors[pore]]  # Neighboring clusters
                unique_ns = np.unique(nc[nc != -1])  # Unique Neighbors
                seq_pore = "S:"+str(un_seq)+" P:"+str(pore)
                if np.all(nc == -1):
                    # This is the start of a new trapped cluster
                    clusters[pore] = next_cluster_num
                    next_cluster_num += 1
                    msg = (seq_pore+" C:1 new cluster number: " +
                           str(clusters[pore]))
                    logger.info(msg)
                elif len(unique_ns) == 1:
                    # Grow the only connected neighboring cluster
                    if not stopped_clusters[unique_ns[0]]:
                        clusters[pore] = unique_ns[0]
                        msg = (seq_pore+" C:2 joins cluster number: " +
                               str(clusters[pore]))
                        logger.info(msg)
                    else:
                        clusters[pore] = -2
                elif -2 in unique_ns:
                    # We have reached a sink neighbor, stop growing cluster
                    msg = (seq_pore+" C:3 joins sink cluster")
                    logger.info(msg)
                    clusters[pore] = -2
                    # Stop growth and merging
                    stopped_clusters[unique_ns[unique_ns > -1]] = True
                else:
                    # We might be able to do some merging
                    # Check if any stopped clusters are neighbors
                    if np.any(stopped_clusters[unique_ns]):
                        msg = (seq_pore+" C:4 joins sink cluster")
                        logger.info(msg)
                        clusters[pore] = -2
                        # Stop growing all neighboring clusters
                        stopped_clusters[unique_ns] = True
                    else:
                        # Merge multiple un-stopped trapped clusters
                        new_num = unique_ns[0]
                        clusters[pore] = new_num
                        for c in unique_ns:
                            clusters[clusters == c] = new_num
                            msg = (seq_pore + " C:5 merge clusters: " +
                                   str(c) + " into "+str(new_num))
                            logger.info(msg)

        # And now return clusters
        self['pore.clusters'] = clusters
        logger.info("Number of trapped clusters" +
                    str(np.sum(np.unique(clusters) >= 0)))
        self['pore.trapped'] = self['pore.clusters'] > -1
        trapped_ts = net.find_neighbor_throats(self['pore.trapped'])
        self['throat.trapped'] = np.zeros([net.Nt], dtype=bool)
        self['throat.trapped'][trapped_ts] = True
        self['pore.invasion_sequence'][self['pore.trapped']] = -1
        self['throat.invasion_sequence'][self['throat.trapped']] = -1
