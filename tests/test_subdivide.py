import OpenPNM
import pytest
print('-----> Using OpenPNM version: '+OpenPNM.__version__)
ctrl = OpenPNM.Base.Controller()
#==============================================================================
'''Build Topological Network'''
#==============================================================================
pn = OpenPNM.Network.Cubic(shape=[5,5,5],spacing=0.001,name='micro_net')
pn['pore.micro'] = True
##============================================================================
nano_pores = [2,13,14,15]
pn.subdivide(pores=nano_pores,shape=[4,4,4],labels='nano')
#============================================================================
assert pn.Nt == (300+(4*144)-16+15*16+16)
##==============================================================================
ctrl.export(network=pn,filename='nano')