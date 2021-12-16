# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 10:28:14 2021

@author: xu kai
"""

import numpy as np

def find_organic_throats(organic_pores, throat_conns, num_throat):
    r"""
    Find organic throats

    Parameters
    ----------
    organic_pores: ndarray
        index of organic pores,
        equuivlant to pores whose diameter is less than 50 nm.

    throat_conns: ndarray, geo['throat.conns']
        pairs of pores which are connected by a throat.

    num_throat: int, geo.Nt
        number of throats

    Return
    ------
    organic_throats: ndarray
        pairs of ndarray
        index 0 is throat conns,
        index 1 is throat labels.

    """

    organic_throats = []
    organic_throat_index = []
    for i in range(num_throat):
        throat = throat_conns[i]
        if (throat[0] in organic_pores and throat[1] in organic_pores):
            organic_throats.append(throat)
            organic_throat_index.append(i)

    organic_throats = np.array(organic_throats)
    organic_throat_index = np.array(organic_throat_index)

    return [organic_throats, organic_throat_index]


