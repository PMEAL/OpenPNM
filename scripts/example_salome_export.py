r"""
In the example script a generic network is created then exported as a
Salome Python script. The script should be executed from Salome with
"load script". The geometry is then built. The geometry generation on
Salome may take some time depending on the number of pores.

"""
import numpy as np
import openpnm as op


# Workspace and project
ws = op.Workspace()
proj = ws.new_project()
export = False

# Network
np.random.seed(7)
net = op.network.Cubic(shape=[4, 3, 3], spacing=1e-4, project=proj)

# Geometry
geo = op.geometry.StickAndBall(network=net, pores=net.Ps, throats=net.Ts)

# Phase
phase = op.phases.Water(network=net)

# Export the network
if export:
    proj.export_data(phases=[phase], filename='out', filetype='Salome')
