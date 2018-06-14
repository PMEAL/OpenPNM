import openpnm as op
from openpnm import topotools
ws = op.Workspace()

bcc = op.network.Bravais(shape=[3, 3, 3], mode='bcc')
fcc = op.network.Bravais(shape=[3, 3, 3], mode='fcc')
fcc['pore.coords'] += [3, 0, 0]

topotools.merge_networks(fcc, bcc)
op.io.VTK.save(network=fcc, filename='sim_02')
