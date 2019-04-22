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
relcalc=op.algorithms.RelativePermeability()
relcalc.setup(input_vect=['xx'])
relcalc.run()