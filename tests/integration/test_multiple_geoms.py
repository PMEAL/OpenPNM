import openpnm as op

pn = op.network.Cubic(shape=[10, 10, 40], spacing=1)

# Generate Ps as boolean mask
Ps = pn['pore.coords'][:, 2] <= 20
Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection', flatten=True)
geom1 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

# Convert Ps to indices
Ps = pn.toindices(pn['pore.coords'][:, 2] > 20)
Ts = pn.find_neighbor_throats(pores=Ps, mode='union')
geom2 = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

water = op.phases.Water(network=pn)

phys1 = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom1)
phys1.add_model(propname='pore.test',
                model=op.models.geometry.pore_misc.constant,
                value=1.0)
phys1.regenerate_models()

phys2 = op.physics.GenericPhysics(network=pn, phase=water, geometry=geom2)
# Copy models from phys1 to phys2 directly
phys2.models = phys1.models
# Ensure the models dicts are equal
assert phys1.models == phys2.models
# But they are actually distinct objects
assert phys1.models is not phys2.models
