import openpnm as op

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4)
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

air = op.phases.Air(network=pn, name='air')
water = op.phases.Water(network=pn, name='h2o')

# phys = op.physics.Standard(network=pn, phase=air, geometry=geo)
# mod = op.models.physics.hydraulic_conductance.stokes_generic
# phys.add_model(propname='throat.hydraulic_conductance_new',
#                model=mod, flow_coeff='throat.flow_coeff')

# sf = op.algorithms.StokesFlow(network=pn, phase=air)
# sf.setup(conductance='throat.hydraulic_conductance')
# sf.set_value_BC(pores=pn.pores('left'), values=101323)
# sf.set_value_BC(pores=pn.pores('right'), values=0)
# sf.run()
# Q1 = sf.rate(pores=pn.pores('left'))
# print(Q1)

# sf = op.algorithms.StokesFlow(network=pn, phase=air)
# sf.setup(conductance='throat.hydraulic_conductance_new')
# sf.set_value_BC(pores=pn.pores('left'), values=101323)
# sf.set_value_BC(pores=pn.pores('right'), values=0)
# sf.run()
# Q2 = sf.rate(pores=pn.pores('left'))
# print(Q2)
