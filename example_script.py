import openpnm as op
from openpnm import topotools
ws = op.Workspace()

bcc = op.network.Bravais(shape=[3, 4, 5], mode='bcc')
fcc = op.network.Bravais(shape=[3, 4, 5], mode='fcc')
fcc['pore.coords'] += [3, 0, 0]
sc = op.network.Bravais(shape=[3, 4, 5], mode='sc')
sc['pore.coords'] += [6, 0, 0]
