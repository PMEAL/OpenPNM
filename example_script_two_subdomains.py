import pytest
import openpnm as op
ws = op.Workspace()
proj = ws.new_project()

pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-4, project=proj)
Ps = pn['pore.coords'][:, 0] < pn['pore.coords'][:, 0].mean()
Ts = pn.find_neighbor_throats(pores=Ps, mode='xnor')
geo1 = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts)

Ps = pn['pore.coords'][:, 0] >= pn['pore.coords'][:, 0].mean()
Ts = pn.find_neighbor_throats(pores=Ps, mode='or')
geo2 = op.geometry.StickAndBall(network=pn, pores=Ps, throats=Ts)

pn['pore.foo'] = 1
# Can't create a subdict below foo
with pytest.raises(Exception):
    pn['pore.foo.bar'] = 1
# Can create a subdict directly
pn['pore.baz.bar'] = 2
# Can't create a new item already used as subdict
with pytest.raises(Exception):
    pn['pore.baz'] = 2

# Also works on subdomains
geo1['pore.blah'] = 1
with pytest.raises(Exception):
    geo1['pore.blah.boo'] = 1
geo1['pore.bee.bop'] = 1
with pytest.raises(Exception):
    geo1['pore.bee'] = 1

# Now start looking across objects
with pytest.raises(Exception):
    geo1['pore.foo'] = 1  # Already exists on pn
with pytest.raises(Exception):
    geo1['pore.foo.bar'] = 1  # pore.foo already exists on pn
with pytest.raises(Exception):
    geo1['pore.baz'] = 1  # pore.baz.bar already exists on pn

# Now start looking across objects
geo2['pore.blah'] = 1
geo2['pore.bee.bop'] = 1
with pytest.raises(Exception):
    geo1['pore.bee'] = 1

with pytest.raises(Exception):
    pn['pore.bee'] = 1

with pytest.raises(Exception):
    pn['pore.bee.bop'] = 1
