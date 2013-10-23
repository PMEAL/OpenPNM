from __future__ import absolute_import
import numpy as np
import OpenPNM

def Cubic(domain_size=[1,1,1], divisions=[10,10,10], stats_pores=None, stats_throats=None, btype=[0,0,0]):
  pn = OpenPNM.Geometry.Cubic().generate(**locals())
  return {'net': pn}
  
def Import(import_file = 'D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles\\large_network'):
    sections = import_file.split('\\')
    filename = sections[np.size(sections,0)-1]
    sections[np.size(sections,0)-1] = ''
    path = sections[0:(np.size(sections,0)-1)].join('\\')
    print( filename )
    print( path )
    pn = OpenPNM.Geometry.MatFile(filename=str(filename),path=str(path)).generate()
    pn.pore_properties['inlets']=np.zeros((pn.get_num_pores(),),dtype=np.int8)
    pn.pore_properties['inlets'][0]=1
    pn.pore_properties['outlets']=np.zeros((pn.get_num_pores(),),dtype=np.int8)
    pn.pore_properties['outlets'][0]=1
    return {'net': pn}