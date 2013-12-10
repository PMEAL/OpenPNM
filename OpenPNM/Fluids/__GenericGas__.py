import OpenPNM
import scipy as sp
import sys
import copy

class GenericGas:
    r"""
    GenericGas - Base class to generate gas properties

    Parameters
    ----------

    """
    def __init__(self,**kwargs):
#        super(GenericFluid,self).__init__(**kwargs)
#        self._logger.debug("Construct class")
        #List of fluid property categories that are invoked when fluid is created
        print('init of GenericGas')

    def create(self,T=298.,P=101325.,**prms):
        r"""
        Create a fluid object using the supplied parameters
        """
        self.props = {}
        key = 'Diffusivity'
        a = getattr(OpenPNM.Fluids, key)
        a = a.__getattribute__(prms[key]['method'])
        self.__setattr__(key,a)
        return self

#        for key in prms:
#            self.__setattr__(key,FluidBuilder(prms[key]))

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
        self.pore_conditions.update({'temperature': 298.})
        self.pore_conditions.update({'pressure': 101325.})
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

        TODO: This method needs plenty of checks,
        - for preexisting pair
        - make sure they're the same temperature and pressure, etc

        '''
        self.partner = fluid2
        fluid2.partner = self

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

    def FluidBuilder():
    
        
        def doit(self,property_name,prms):
            self.params = {}
            for it in prms.keys():
                self.params[it] = prms[it]
            a = getattr(OpenPNM.Fluids, property_name)
            self.method = a.__getattribute__(prms['method'])
            return self

if __name__ =="__main__":

    #Create fluids
    agg = OpenPNM.Fluids.GenericFluid().create(agg_recipe)


