import openpnm as op
from openpnm import topotools
ws = op.Workspace()

bcc = op.network.Bravais(shape=[4, 5, 6], mode='bcc')
fcc = op.network.Bravais(shape=[4, 5, 6], mode='fcc')
fcc['pore.coords'] += [4, 0, 0]

topotools.merge(fcc, bcc)
op.io.VTK.save(network=fcc, filename='sim_02')
