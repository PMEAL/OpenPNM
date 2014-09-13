import OpenPNM
import OpenPNM.Geometry.models as gm

class Sandstone():
    
    @staticmethod
    def Berea(shape,name=''):
        r'''
        '''
        name = 'Berea'+name
        #Generate Cubic topology with correct lattice spacing
        pn = OpenPNM.Network.Cubic(shape=shape,name=name,spacing=0.000040)
        
        #Generate a GenericGeometry object
        geom = OpenPNM.Geometry.GenericGeometry(network=pn)
        
        #Add appropriate models to Geometry
        geom.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       num_range=[0,0.95],
                       seed=None)
        
        return pn
        
    @staticmethod
    def Fountainebleu(shape,name=None):
        r'''
        '''
        name = 'SGL10BA'+name
        #Generate Cubic topology with correct lattice spacing
        pn = OpenPNM.Network.Cubic(shape=shape,name=name,spacing=0.000040)
        
        #Generate a GenericGeometry object
        geom = OpenPNM.Geometry.GenericGeometry(network=pn)
        
        #Add appropriate models to Geometry
        geom.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       num_range=[0,0.8834],
                       seed=None)
        
        return pn