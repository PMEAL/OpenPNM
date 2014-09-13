import OpenPNM
import OpenPNM.Geometry.models as gm

class GDL():
    
    @staticmethod
    def Toray090_Cubic(shape,name=''):
        r'''
        '''
        name = 'Toray090'+name
        #Generate Cubic topology with correct lattice spacing
        pn = OpenPNM.Network.Cubic(shape=shape,name=name,spacing=0.000040)
        
        #Generate a GenericGeometry object
        geom = OpenPNM.Geometry.GenericGeometry(network=pn)
        
        #Add appropriate models to Geometry
        geom.add_model(propname='pore.seed',
                       model=gm.pore_misc.random,
                       num_range=[0,0.95],
                       seed=None)
        geom.remove_model('pore.seed')
        geom.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        geom.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.sphere,
                       psd_name='weibull_min',
                       psd_shape=2.77,
                       psd_loc=6.9e-7,
                       psd_scale=9.8e-6,
                       psd_offset=10e-6)
        geom.add_model(propname='pore.area',
                       model=gm.pore_area.spherical)
        geom.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        geom.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.cylinder,
                       tsd_name='weibull_min',
                       tsd_shape=2.77,
                       tsd_loc=6.9e-7,
                       tsd_scale=9.8e-6,
                       tsd_offset=10e-6)             
        geom.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        geom.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)
        geom.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder)
        geom.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)
        
        return pn
        
    @staticmethod
    def SGL10BA_Cubic(shape,name=None):
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
        geom.remove_model('pore.seed')
        geom.add_model(propname='throat.seed',
                       model=gm.throat_misc.neighbor,
                       pore_prop='pore.seed',
                       mode='min')
        geom.add_model(propname='pore.diameter',
                       model=gm.pore_diameter.sphere,
                       psd_name='weibull_min',
                       psd_shape=3.07,
                       psd_loc=1.97e-6,
                       psd_scale=1.6e-5,
                       psd_offset=18e-6)
        geom.add_model(propname='pore.area',
                       model=gm.pore_area.spherical)
        geom.add_model(propname='pore.volume',
                       model=gm.pore_volume.sphere)
        geom.add_model(propname='throat.diameter',
                       model=gm.throat_diameter.cylinder,
                       tsd_name='weibull_min',
                       tsd_shape=3.07,
                       tsd_loc=1.97e-6,
                       tsd_scale=1.6e-5,
                       tsd_offset=18e-6)               
        geom.add_model(propname='throat.length',
                       model=gm.throat_length.straight)
        geom.add_model(propname='throat.volume',
                       model=gm.throat_volume.cylinder)
        geom.add_model(propname='throat.area',
                       model=gm.throat_area.cylinder)
        geom.add_model(propname='throat.surface_area',
                       model=gm.throat_surface_area.cylinder)
        
        return pn