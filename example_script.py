import openpnm as op
ws = op.Workspace()
proj = ws.new_project()

pn = op.network.Cubic(shape=[15, 15, 15], spacing=1e-4, project=proj)
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phases.Air(network=pn)
water = op.phases.Water(network=pn)
hg = op.phases.Mercury(network=pn)
phys_air = op.physics.Standard(network=pn, phase=air, geometry=geo)
phys_water = op.physics.Standard(network=pn, phase=water, geometry=geo)
phys_hg = op.physics.Standard(network=pn, phase=hg, geometry=geo)

mip = op.algorithms.Porosimetry(network=pn)
mip.setup(phase=hg)
mip.set_inlets(pores=pn.pores(['top', 'bottom']))
mip.run()
# mip.plot_intrusion_curve()

perm = op.algorithms.StokesFlow(network=pn)
perm.setup(phase=water)
perm.set_value_BC(pores=pn.pores('right'), values=0)
perm.set_value_BC(pores=pn.pores('left'), values=101325)
perm.run()
# print(perm.calc_eff_permeability())

# Add reaction term to phys_air
phys_air['pore.n'] = 2
phys_air['pore.A'] = -1e-5
phys_air.add_model(propname='pore.2nd_order_rxn',
                   model=op.models.physics.generic_source_term.standard_kinetics,
                   quantity='pore.concentration',
                   prefactor='pore.A', exponent='pore.n')
rxn = op.algorithms.FickianDiffusion(network=pn)
rxn.setup(phase=air)
Ps = pn.find_nearby_pores(pores=2000, r=5e4)
rxn.set_source(propname='pore.2nd_order_rxn', pores=Ps)
rxn.set_value_BC(pores=pn.pores('top'), values=1)
rxn.run()
