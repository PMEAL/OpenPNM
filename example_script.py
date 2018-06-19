import openpnm as op
ws = op.Workspace()
proj = ws.new_project()
pn = op.network.Cubic(shape=[5, 5, 5], project=proj)
Ps = pn.pores(['left', 'right'])
Ps = pn.tomask(Ps)
geom1 = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phases.Air(network=pn)
phys1 = op.physics.Standard(network=pn, phase=air, geometry=geom1)

alg = op.algorithms.GenericTransport(network=pn)
alg.setup(phase=air, quantity='pore.mole_fraction',
          conductance='throat.diffusive_conductance')
alg.set_value_BC(pores=pn.pores('left'), values=1)
alg.set_value_BC(pores=pn.pores('right'), values=0)
alg.run()
