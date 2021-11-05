# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 09:08:50 2021

@author: xu kai

from https://docs.scipy.org/doc/scipy/reference/tutorial/stats.html
"""

import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from scipy import stats

def bimodal_distribution(pore_num):
    r"""
    generate bimodal distribution for pore diameters
    in this function, I set the pore size diameter in the range of 25 - 250 microns
    """

    np.random.seed(1)
    num1 = int(pore_num * 0.2)
    num2 = pore_num - num1
    loc1, scale1, size1 = (-2, 0.2, num1)
    loc2, scale2, size2 = (2, 1, num2)

    x2 = np.concatenate([np.random.normal(loc=loc1, scale=scale1, size=size1),
                         np.random.normal(loc=loc2, scale=scale2, size=size2)])

#    x_eval = np.linspace(x2.min() - 1, x2.max() + 1, 1000)
#
#    def my_kde_bandwidth(obj, fac=1./5):
#        """We use Scott's Rule, multiplied by a constant factor."""
#        return np.power(obj.n, -1./(obj.d+4)) * fac
#
#    kde = stats.gaussian_kde(x2)
#    kde2 = stats.gaussian_kde(x2, bw_method='silverman')
#    kde3 = stats.gaussian_kde(x2, bw_method=partial(my_kde_bandwidth, fac=0.2))
#    kde4 = stats.gaussian_kde(x2, bw_method=partial(my_kde_bandwidth, fac=0.5))
#
#    pdf = stats.norm.pdf
#    bimodal_pdf = pdf(x_eval, loc=loc1, scale=scale1) * float(size1) / x2.size + \
#                  pdf(x_eval, loc=loc2, scale=scale2) * float(size2) / x2.size
#
#    fig = plt.figure(figsize=(8, 6))
#    ax = fig.add_subplot(111)
#
#    ax.plot(x2, np.zeros(x2.shape), 'b+', ms=12)
#    ax.plot(x_eval, kde(x_eval), 'k-', label="Scott's Rule")
#    ax.plot(x_eval, kde2(x_eval), 'b-', label="Silverman's Rule")
#    ax.plot(x_eval, kde3(x_eval), 'g-', label="Scott * 0.2")
#    ax.plot(x_eval, kde4(x_eval), 'c-', label="Scott * 0.5")
#    ax.plot(x_eval, bimodal_pdf, 'r--', label="Actual PDF")
#
#    ax.set_xlim([x_eval.min(), x_eval.max()])
#    ax.legend(loc=2)
#    ax.set_xlabel('x')
#    ax.set_ylabel('Density')
#    plt.show()

    # rescale pore size to 25 micros to 250 micros
    x2 = (x2 - x2.min()) / (x2.max() - x2.min()) * (250e-9 - 25e-9) + 25e-9
    np.random.shuffle(x2)

    return x2

