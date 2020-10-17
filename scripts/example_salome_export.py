import openpnm as op
import numpy as np

# work space and project
ws = op.Workspace()
ws.settings["loglevel"] = 30
proj = ws.new_project()

# network
np.random.seed(7)
net = op.network.Cubic(shape=[4, 3, 3], spacing=1e-4, project=proj)

# geometry
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)

# phase
phase = op.phases.Water(network=net)

# export the network
proj.export_data(phases=[phase], filename='OUT', filetype='Salome')
