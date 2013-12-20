
"""
module Physics
===============================================================================

"""
import OpenPNM
import scipy as sp
from functools import partial

class GenericFluid(OpenPNM.Base.Utilities):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        self.pore_conditions = {}
        self.throat_conditions = {}

    def create(self,network,T=298.,P=101325.,**recipe):
        r"""
        Create a fluid object using the supplied parameters
        """
        try: recipe = self.recipe #check if recipe is pre-existing on self (from init of subclassed methods)
        except: pass
        try: self.name = recipe['name']
        except: self._logger.error('Fluid name must be given')
        network._fluids.append(self)
        self.Tc = recipe['Tc']
        self.Pc = recipe['Pc']
        self.MW = recipe['MW']
        self.pore_conditions['temperature'] = T
        self.pore_conditions['pressure'] = P
        for key, args in recipe.items():
            try:
                function = getattr( getattr(OpenPNM.Fluids, key), args['method'] ) #Get method from the file
                preloaded_fn = partial(function, fluid=self, network=network, **args) 
                setattr(self, key, preloaded_fn)
                self._logger.info('Successfully added '+key+' to '+self.name)
            except AttributeError: pass
        self.regenerate()
        return self

    def regenerate(self):
        r'''
        This updates all properties using the methods indicated in the recipe.
        '''
        try: self.viscosity()
        except: pass
        try: self.diffusivity()
        except: pass
        try: self.molar_density()
        except: pass
        try: self.surface_tension()
        except: pass
        try: self.contact_angle()
        except: pass

    def interpolate_pore_conditions(self,network,Tinfo=None):
        r"""
        Determines a pore property as the average of it's neighboring throats

        Parameters
        ----------
        network : OpenPNM Pore Network Object
            The network on which to perform the interpolation

        Tinfo : array_like
            The array of throat information to be interpolated

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.

        """
        if sp.size(Tinfo)==1:
            Pinfo = Tinfo
        elif sp.size(Tinfo) != network.get_num_throats():
            raise Exception('The list of throat information received was the wrong length')
        else:
            Pinfo = sp.zeros((network.get_num_pores()))
            #Only interpolate conditions for internal pores, type=0
            Pnums = sp.r_[0:network.get_num_pores(Ptype=[0])]
            nTs = network.get_neighbor_throats(Pnums,flatten=False)
            for i in sp.r_[0:sp.shape(nTs)[0]]:
                Pinfo[i] = sp.mean(Tinfo[nTs[i]])
        return Pinfo

    def interpolate_throat_conditions(self,network,Pinfo=None):
        r"""
        Determines a throat condition as the average of the conditions it's neighboring pores

        Parameters
        ----------
        network : OpenPNM Pore Network Object
            The network on which to perform the interpolation

        Pinfo : array_like
            The name of the throat condition to be interpolated

        Notes
        -----
        This uses an unweighted average, without attempting to account for distances or sizes of pores and throats.

        """
        if sp.size(Pinfo)==1:
            Tinfo = Pinfo
        elif sp.size(Pinfo) != network.get_num_pores():
            raise Exception('The list of pore information received was the wrong length')
        else:
            Tinfo = sp.zeros((network.get_num_throats()))
            #Interpolate values for all throats, including those leading to boundary pores
            Tnums = sp.r_[0:network.get_num_throats()]
            nPs = network.get_connected_pores(Tnums,flatten=False)
            for i in sp.r_[0:sp.shape(nPs)[0]]:
                Tinfo[i] = sp.mean(Pinfo[nPs[i]])
        return Tinfo
        
    def __str__(self):
        return('test')

if __name__ =="__main__":

    #Create fluids
    air_recipe = {
    'name': 'air',
    'Pc': 3.771e6, #Pa
    'Tc': 132.65,  #K
    'MW': 0.0291,  #kg/mol
    'diffusivity': {'method': 'Fuller',
                    'MA': 0.03199,
                    'MB': 0.0291,
                    'vA': 16.3,
                    'vB': 19.7},
    'viscosity': {'method': 'Reynolds',
                  'uo': 0.001,
                  'b': 0.1},
    'molar_density': {'method': 'ideal_gas',
                      'R': 8.314},
    'surface_tension': {'method': 'constant',
                        'value': 0},
    'contact_angle': {'method': 'na'},
    }
    gas = OpenPNM.Fluids.GenericGas().create(air_recipe)



