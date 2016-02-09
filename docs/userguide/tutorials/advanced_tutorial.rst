


>>> geom1.models.add(propname='pore.seed',
...                  model=OpenPNM.Geometry.models.pore_seed.random,
...                  seed=1,
...                  num_range=[0, 1],
...                  regen_mode='constant')
>>> geom2.models['pore.seed'] = geom1.models['pore.seed'].copy()
