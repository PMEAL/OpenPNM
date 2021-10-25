import openpnm as op
from openpnm.materials import BundleOfTubes


ws = op.Workspace()
ws.settings['loglevel'] = 30
proj = ws.new_project()

net, geo, phase = BundleOfTubes(
    shape=50,
    spacing=0.0001,
    length=0.001,
    psd_params={
        'distribution': 'weibull',
        'loc': 1e-6,
        'scale': 4e-5,
        'shape': 2.2
    },
    settings={'adjust_psd': 'normalize', 'seed': 0}
)

hg = op.phases.Mercury(network=net)
phys = op.physics.Classic(network=net, phase=hg, geometry=geo)

mip = op.algorithms.Porosimetry(network=net, phase=hg)
mip.set_inlets(pores=net.pores('top'))
mip.run()
mip.plot_intrusion_curve()
