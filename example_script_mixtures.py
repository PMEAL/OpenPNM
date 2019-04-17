import openpnm as op
ws = op.Workspace()
proj = ws.new_project()

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4, project=proj)
geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)

N2 = op.phases.components.gases.N2(network=pn, name='pure_N2')
O2 = op.phases.components.gases.O2(network=pn, name='pure_O2')
air = op.phases.GenericMixture(network=pn, components=[N2, O2],
                               name='air_mixture')
air.set_mole_fraction(N2, 0.791)
air.set_mole_fraction(O2, 0.209)
air.add_model(propname='pore.molar_mass',
              model=op.models.phases.mixtures.mole_weighted_average,
              prop='pore.molecular_weight')
air.add_model(propname='pore.diffusivity',
              model=op.models.phases.mixtures.fuller)
air.add_model(propname='pore.viscosity',
              model=op.models.misc.polynomial,
              prop='pore.temperature',
              a=[0.00000182082, 6.51815E-08, -3.48553E-11,
                 1.11409E-14])

phys = op.physics.GenericPhysics(network=pn, phase=air, geometry=geo)
phys.add_model(propname='throat.diffusive_conductance',
               model=op.models.physics.diffusive_conductance.ordinary_diffusion)

fd = op.algorithms.FickianDiffusion(network=pn, phase=air)
fd.set_value_BC(pores=pn.pores('left'), values=1)
fd.set_value_BC(pores=pn.pores('right'), values=0)
fd.run()
