import numpy as np
import openpnm as op

pn = op.core.FullDomain()

# The number of allowed pores and throats is dictated by the initial arrays
pn['pore.anything'] = np.ones(6)
pn['throat.anything'] = np.ones(6)
print(pn.count('pore'))
print(pn.count('throat'))

# Let's make two geometry objects
geo1 = op.core.SubDomain(domain=pn, name='geo1', pores=[0, 1, 4], throats=[0, 1, 5])
geo2 = op.core.SubDomain(domain=pn, name='geo2', pores=[2, 3, 5], throats=[2, 3, 4])

# Now let's inspect each object
print(pn)
print(geo1)  # This is just a view into the actual array on pn
print(geo2)

# Let's write data to one geo and watch interleave data
geo1['pore.new'] = [11, 22, 33]
print(geo1)
print(geo2['pore.new'])
print(pn['pore.new'])

# We can also assign scalars and have them broadcast
geo2['pore.scalar'] = 2
# And like old behavior, the empty spots are nans, so array is float
print(pn['pore.scalar'])

# Let's move pore 4 onto geo2
pn['pore.geo1'][4] = False
pn['pore.geo2'][4] = True
# Not only is the value missing from geo1, it's now part of geo2
print(geo1['pore.new'])
print(geo2['pore.new'])

# Deleting pores and/or throats is super easy!
pn._data['pore'].drop([4])
print(geo1['pore.new'])
print(geo2['pore.new'])

# We can ask each subdomain where it applies
print(geo1.locs('pore'))
print(geo2.locs('pore'))

# A 'param' prefix is now supported
pn['param.k'] = 1.4
print(pn)
# Note that user specified prefixes won't work since we won't know how to
# deal with them when deleting pores for instance, so we will stick with
# pore, throat, and prefix, and maybe conduit?

# Labels still work on SubDomains
geo2['pore.label'] = [True, True, False]
print(geo2.pores('label'))
# But are also available on the other domain, and full domain
print(geo1.pores('label'))
print(pn.pores('label'))
# The interesting part is the new 'relative' keyword
print(geo2.pores('label', relative=False))  # No more mapping!

# The following does NOT work any more:
print(geo1['pore.new'])
geo1['pore.new'][0] = 55
print(geo1['pore.new'])
# For this we have the set_values function
geo1.set_values(key='pore.new', indices=[0], values=[55])
print(geo1['pore.new'])
# Or just use the full domain object
pn['pore.new'][0] = 66
print(geo1['pore.new'])
