import matplotlib.pyplot as plt
import openpnm as op
import scipy as sp
ws = op.Workspace()
proj = ws.new_project()

pn = op.network.Bravais(shape=[6, 5, 4], mode='fcc', project=proj)

op.io.VTK.save(network=pn, filename='sim_02')
