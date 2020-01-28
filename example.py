import openpnm as op
ws = op.Workspace()
proj = ws.new_project()

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4, project=proj)
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
air = op.phases.Air(network=pn, name='air')
water = op.phases.Water(network=pn, name='h2o')
hg = op.phases.Mercury(network=pn, name='hg')
phys_air = op.physics.Standard(network=pn, phase=air, geometry=geo)
phys_water = op.physics.Standard(network=pn, phase=water, geometry=geo)
phys_hg = op.physics.Standard(network=pn, phase=hg, geometry=geo)

mip = op.algorithms.Porosimetry(network=pn)
mip.setup(phase=hg)
mip.set_inlets(pores=pn.pores(['top', 'bottom']))
mip.run()
hg.update(mip.results(Pc=70000))
# mip.plot_intrusion_curve()

perm = op.algorithms.StokesFlow(network=pn)
perm.setup(phase=water)
perm.set_value_BC(pores=pn.pores('right'), values=0)
perm.set_value_BC(pores=pn.pores('left'), values=101325)
perm.run()
water.update(perm.results())
# print(perm.calc_effective_permeability())

# Add reaction term to phys_air
mod = op.models.physics.generic_source_term.standard_kinetics
phys_air['pore.n'] = 2
phys_air['pore.A'] = -1e-5
phys_air.add_model(propname='pore.2nd_order_rxn', model=mod,
                   quantity='pore.concentration',
                   prefactor='pore.A', exponent='pore.n',
                   regen_mode='deferred')
rxn = op.algorithms.FickianDiffusion(network=pn)
rxn.setup(phase=air, solver='spsolve')
Ps = pn.find_nearby_pores(pores=50, r=5e-4, flatten=True)
rxn.set_source(propname='pore.2nd_order_rxn', pores=Ps)
rxn.set_value_BC(pores=pn.pores('top'), values=1)
rxn.run()
air.update(rxn.results())

fd = op.algorithms.FickianDiffusion(network=pn)
fd.setup(phase=air)
fd.set_value_BC(pores=pn.pores('left'), values=1)
fd.set_value_BC(pores=pn.pores('right'), values=0)
fd.run()
fd.calc_effective_diffusivity()

# Output network and phases to a VTP file for visualization in Paraview
# proj.export_data(network=pn, phases=[hg, air, water], filename='output.vtp')
