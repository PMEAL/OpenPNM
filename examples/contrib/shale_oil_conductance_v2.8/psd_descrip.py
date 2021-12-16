# -*- coding: utf-8 -*-
"""
Created on Sat Jul 10 20:27:46 2021

@author: xu kai
"""

import matplotlib.pyplot as plt

def psd_descrip(pn, x, y, z):
    """
    pnm description, including
    * sample size in x y z directions;
    * porosity;
    * pore size and throat size;

    x y z are axis lengths.
    """
    # sample size
    pn['pore.coords'][:,0].min()
    pn['pore.coords'][:,0].max()

    pn['pore.coords'][:,1].min()
    pn['pore.coords'][:,1].max()

    pn['pore.coords'][:,2].min()
    pn['pore.coords'][:,2].max()

    len_x = pn['pore.coords'][:, 0].max() - pn['pore.coords'][:, 0].min()
    len_y = pn['pore.coords'][:, 1].max() - pn['pore.coords'][:, 1].min()
    len_z = pn['pore.coords'][:, 2].max() - pn['pore.coords'][:, 2].min()

    pn['pore.diameter'] = 2 * pn['pore.radius']
    pn['throat.diameter'] = 2 * pn['throat.radius']

    # pore size
    pd_max = pn['pore.diameter'].max()
    pd_mean = pn['pore.diameter'].mean()
    pd_min = pn['pore.diameter'].min()

    # throat size
    td_max = pn['throat.diameter'].max()
    td_mean = pn['throat.diameter'].mean()
    td_min = pn['throat.diameter'].min()

    # porosity calculation
    v = (x * y * z)
    p1 = (pn['pore.volume'].sum() + pn['throat.volume'].sum()) / v
    p2 = (pn['pore.volume'].sum()) / v

    psd_dict = {'sample size': (len_x, len_y, len_z),
                'porosity with throats': p1,
                'porosity no throats': p2,
                'pore dia max': pd_max,
                'pore dia mean': pd_mean,
                'pore dia min': pd_min,
                'throat dia max': td_max,
                'throat dia mean': td_mean,
                'throat dia min': td_min}

    print(psd_dict)

    # plot
    pore_dia=pn['pore.diameter'],
    throat_dia=pn['throat.diameter'],
    throat_len=pn['throat.length'],

    fig2, (ax1, ax2) = plt.subplots(1, 2)

    # pore_diameter, throat_diameter
    ax1.hist(pore_dia, density=False, histtype='stepfilled', bins=50,
             color='b', alpha=0.2, label='pore diameter ', edgecolor='k')
    ax1.hist(throat_dia, density=False, bins=50, histtype='stepfilled',
             color='r', alpha=0.2, label='throat diameter', edgecolor='k')

    # throat_length
    ax2.hist(throat_len, density=False,
             histtype='stepfilled', color='g', alpha=0.2, label='throat length',
             edgecolor='k', bins=50)

    ax1.set_title('pore and throat diameter diatribution')
    ax2.set_title('throat length diatribution')



