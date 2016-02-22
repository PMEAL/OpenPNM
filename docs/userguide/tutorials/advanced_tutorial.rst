.. _advanced_tutorial:

###############################################################################
Tutorial 3 of 3: Advaned Topics and Usage
###############################################################################

**Learning Outcomes**

1. Employ data exchange bewteen objects
2. Combine multiple algorithms
3. Define Physics with pores and throats instead of Geometry
4. Explore the ModelsDict design
5. Use the workspace manager to save and load


>>> geom1.models.add(propname='pore.seed',
...                  model=OpenPNM.Geometry.models.pore_seed.random,
...                  seed=1,
...                  num_range=[0, 1],
...                  regen_mode='constant')
>>> geom2.models['pore.seed'] = geom1.models['pore.seed'].copy()
