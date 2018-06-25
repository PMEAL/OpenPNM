import openpnm as op
ws = op.Workspace()
proj = ws.new_project()
pn = op.network.Cubic(shape=[5, 5, 5], project=proj)
geom = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phases.Air(network=pn)
phys = op.physics.Standard(network=pn, phase=air, geometry=geom)

alg = op.algorithms.GenericTransport(network=pn)
alg.setup(phase=air, quantity='pore.mole_fraction',
          conductance='throat.diffusive_conductance')
alg.set_value_BC(pores=pn.pores('left'), values=1)
alg.set_value_BC(pores=pn.pores('right'), values=0)
alg.run()

for item in air.__dir__():
    if not item.startswith('_'):
        if item not in dict().__dir__():
            if item not in op.core.Base().__dir__():
                print(item)
