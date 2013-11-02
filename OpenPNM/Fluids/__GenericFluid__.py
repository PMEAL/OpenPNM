import OpenPNM
import scipy as sp
import copy

class GenericFluid(OpenPNM.Utilities.OpenPNMbase):
    r"""
    GenericFluid - Base class to generate fluid properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
        super(GenericFluid,self).__init__(**kwargs)
        self._logger.debug("Construct class")
        #List of fluid property categories that are invoked when fluid is created
        self._implemented_methods = ['diffusivity',
                                     'viscosity',
                                     'molar_density',
                                     'surface_tension',
                                     'contact_angle']

    def create(self,fluid_recipe,T=298.,P=101325.):
        r"""
        Creates a fluid using the recipe provided

        Parameters
        ----------
        fluid_recipe : Dictionary
            A dictionary of key : value pairs that provide instructions for calculating fluid conditions

        T & P : float (optional)
            The temperature and pressure at which fluid should be create, defaults to STP.
        """
        self._fluid_recipe = copy.deepcopy(fluid_recipe)
        self.pore_conditions = {}
        self.throat_conditions = {}
        self.pore_conditions.update({'temperature': T})
        self.pore_conditions.update({'pressure': P})
        self.regenerate()
        return self

    def refresh(self):
        self.regenerate()

    def regenerate(self):
        r'''
        This updates all properties using the methods indicated in the recipe.  This method also takes the opportunity to ensure all values are Numpy arrays.
        '''
        for condition in self._implemented_methods:
            self.pore_conditions.update({condition: getattr(self,condition)()})
        #Make sure all values are Numpy arrays (Np,) or (Nt,)
        for i in self.pore_conditions.keys():
            self.pore_conditions[i] = sp.array(self.pore_conditions[i],ndmin=1)
        for i in self.throat_conditions.keys():
            self.throat_conditions[i] = sp.array(self.throat_conditions[i],ndmin=1)

    def reset(self):
        r'''
        Remove all existing condtions from the fluid

        TODO: This works, but is kludgy
        '''
        self.pore_conditions = {}
        self.throat_conditions = {}
        try: del self.partner
        except: pass
        self.pore_conditions.update({'temperature': 298})
        self.pore_conditions.update({'pressure': 101325})
        self.regenerate()
        return self

    def clone(self):
        r'''
        Create an exact duplicate fluid, but a unique object.

        TODO: Doesn't work yet
        '''
        return self.__new__

    def set_pair(self,fluid2):
        r'''
        Creates a fluid pair by storing each fluid object on the other.  This allows tracking of defending vs invading phase, among other things.

        TODO: This method needs plenty of checks, for preexisting pair, etc.

        '''
        self.partner = fluid2
        fluid2.partner = self

    def diffusivity(self):
        params = self._fluid_recipe['diffusivity']
        eqn = getattr(OpenPNM.Fluids.Diffusivity,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

    def viscosity(self):
        params = self._fluid_recipe['viscosity']
        eqn = getattr(OpenPNM.Fluids.Viscosity,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

    def molar_density(self):
        params = self._fluid_recipe['molar_density']
        eqn = getattr(OpenPNM.Fluids.MolarDensity,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

    def surface_tension(self):
        params = self._fluid_recipe['surface_tension']
        eqn = getattr(OpenPNM.Fluids.SurfaceTension,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

    def contact_angle(self):
        params = self._fluid_recipe['contact_angle']
        eqn = getattr(OpenPNM.Fluids.ContactAngle,params['method'])
        vals = eqn(self,**params)
        return sp.array(vals,ndmin=1)

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

if __name__ =="__main__":

    #Create fluids
    air = OpenPNM.Fluids.Air().create(T=353,P=200000)

    print ''
    print 'current pore conditions:'
    for i in air.pore_conditions.keys():
        print i,'=',air.pore_conditions[i]
    print ''
