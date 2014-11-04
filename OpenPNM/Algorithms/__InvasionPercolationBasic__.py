# -*- coding: utf-8 -*-
"""
===============================================================================
InvasionPercolationBasic: Simple IP
===============================================================================

"""

import scipy as sp
import bisect
from OpenPNM.Algorithms import GenericAlgorithm

class InvasionPercolationBasic(GenericAlgorithm):
    r"""
    Blah

    Note
    ----
    n/a

    """

    def __init__(self,**kwargs):
        r'''

        '''
        super(InvasionPercolationBasic,self).__init__(**kwargs)

    def run(self,phase,inlets):
        r'''
        '''
        phase = phase
        net = self._net
        # Setup arrays and info
        t_entry = phase['throat.capillary_pressure']
        t_sorted = sp.argsort(t_entry,axis=0)  # Indices into t_entry giving a sorted list
        t_order = sp.zeros_like(t_sorted)
        t_order[t_sorted] = sp.arange(0,net.Nt)  # Location in sorted list
        t_inv = -sp.ones_like(net.Ts)  # List for tracking throat invasion order
        t_seq = -sp.ones_like(net.Ts)  # List for recreating invasion sequence
        p_inv = -sp.ones_like(net.Ps)  # List for tracking pore invasion order
        p_seq = -sp.ones_like(net.Ps)  # List for recreating invasion sequence
        p_inv[inlets] = 0  # Set inlet pores to invaded
        # Perform initial analysis on input pores
        Ts = net.find_neighbor_throats(pores=inlets)
        queue = []
        [bisect.insort(queue,t) for t in t_order[Ts]]  # Using binary search tree for sorting
        tcount = 0
        pcount = 0
        while len(queue) > 0:
            t = queue.pop(0)  # Find throat at the top of the queue
            t_next = t_sorted[t]  # Extract actual throat number
            t_inv[t_next] = tcount  # Note invasion sequence
            t_seq[tcount] = t_next
            while (len(queue)>0) and (queue[0] == t):  # Ensure throat is not duplicated
                t = queue.pop(0)  # Note: Preventing duplicate entries below might save some time here
            Ps = net['throat.conns'][t_next]  # Find pores connected to newly invaded throat
            Ps = Ps[p_inv[Ps]<0]  # Remove already invaded pores from Ps
            if len(Ps)>0:
                pcount += 1
                p_inv[Ps] = pcount  # Note invasion sequence
                p_seq[pcount] = Ps
                Ts = net.find_neighbor_throats(pores=Ps)  # Find connected throats
                Ts = Ts[t_inv[Ts]<0]  # Remove already invaded throats from Ts
                [bisect.insort(queue,t) for t in t_order[Ts]]  # Add new throats to queue
            tcount += 1





if __name__ == '__main__':
    print('no tests yet')
