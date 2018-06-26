import openpnm as op
from openpnm import topotools
ws = op.Workspace()

bcc = op.network.Bravais(shape=[3, 4, 5], mode='bcc')
fcc = op.network.Bravais(shape=[3, 4, 5], mode='fcc')
fcc['pore.coords'] += [3, 0, 0]
sc = op.network.Bravais(shape=[3, 4, 5], mode='sc')
sc['pore.coords'] += [6, 0, 0]

for item in air.__dir__():
    if not item.startswith('_'):
        if item not in dict().__dir__():
            if item not in op.core.Base().__dir__():
                print(item)
