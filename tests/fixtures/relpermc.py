# -*- coding: utf-8 -*-
"""
Created on Wed Jan  9 14:59:50 2019

@author: work
"""
# sample code for relperm
import openpnm as op
import matplotlib.pyplot as plt
import numpy as np
pn = op.network.Cubic(shape=[20, 20, 20], spacing=0.00006)
geom = op.geometry.StickAndBall(network=pn, pores=pn['pore.all'],
                                throats=pn['throat.all'])
relcalc=op.algorithms.RelativePermeability(network=pn)
relcalc.setup(input_vect=['xx','xy','xz'])
relcalc.plot_rel_perm()
#relcalc.setup_ip_algs()
#named_tuple_of_data = relcalc.get_krel_data()
#fig = relcalc.plot_krel_curves()
#named_tuple_of_fitted_exponents = relcalc.fit_krel_curves()
