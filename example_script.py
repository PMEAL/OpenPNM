import matplotlib.pyplot as plt
import openpnm as op
import scipy as sp
ws = op.Workspace()

bcc = op.network.Bravais(shape=[2, 3, 4], mode='bcc')
fcc = op.network.Bravais(shape=[2, 3, 4], mode='fcc')

op.io.VTK.save(network=fcc, filename='sim_02')
