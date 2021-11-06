# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 20:11:18 2021

@author: xu kai
"""

import numpy as np
import matplotlib.pyplot as plt

def pore_coordination_num(pn):
    '''
    calculate pore coordination number
    parameters:
        throats=pn['throat.conns']: array,(pn.Nt, 2), throat with their connected pores;
        pore_num=pn.Np: nubmer of pores
    '''
    throats=pn['throat.conns']
    pore_num=pn.Np
    coordination = np.zeros(pore_num)

    for throat in throats:
        coordination[throat[0]] += 1
        coordination[throat[1]] += 1

    coordination = np.array(coordination, dtype=int)
    bins = np.max(coordination) + 1
    bin_range = (coordination.min()-0.5, coordination.max()+0.5)

    plt.hist(coordination, bins=bins, align='mid', range=bin_range,
             edgecolor='black')
    plt.show()
    return (coordination)  # type change to int number
