"""
module __GenericGeometry__: Base class to construct pore networks
==================================================================

"""

import OpenPNM
from OpenPNM.Base import Core
import scipy as sp

class GenericGeometry(Core):
    r"""
    GenericGeometry - Base class to construct a Geometry object

    Parameters
    ----------
    network : OpenPNM Network Object

    name : string
        A unique name to apply to the object.  This name will also be used as a
        label to identify where this this geometry applies.

    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)

    loggername : string (optional)
        Sets a custom name for the logger, to help identify logger messages

    Examples
    --------
    >>> pn = OpenPNM.Network.TestNet()
    >>> Ps = pn.pores()  # Get all pores
    >>> Ts = pn.throats()  # Get all throats
    >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=Ps,throats=Ts)
    """

    def __init__(self,network=None,pores=[],throats=[],name=None,**kwargs):
        r"""
        Initialize
        """
        super(GenericGeometry,self).__init__(**kwargs)
        self._logger.debug("Class Constructor")

        #Initialize locations
        self['pore.all'] = sp.array([],ndmin=1,dtype=bool)
        self['throat.all'] = sp.array([],ndmin=1,dtype=bool)

        if network == None:
            self._net = OpenPNM.Network.GenericNetwork()
        else:
            self._net = network  # Attach network to self
            self._net._geometries.append(self)  # Register self with network.geometries
        self.name = name

        #Initialize a label dictionary in the associated network
        self._net['pore.'+self.name] = False
        self._net['throat.'+self.name] = False
        self.set_locations(pores=pores,throats=throats)

    def set_locations(self,pores=[],throats=[]):
        r'''
        '''
        if pores != []:
            #Initialize locations
            self['pore.all'] = sp.ones((sp.shape(pores)[0],),dtype=bool)
            self['pore.map'] = pores
            #Specify Geometry locations in Network dictionary
            self._net['pore.'+self.name][pores] = True
        if throats != []:
            #Initialize locations
            self['throat.all'] = sp.ones((sp.shape(throats)[0],),dtype=bool)
            self['throat.map'] = throats
            #Specify Geometry locations in Network dictionary
            self._net['throat.'+self.name][throats] = True

    def check_geometry_health(self):
        r'''
        Perform a check to find pores with overlapping or undefined Geometries
        '''
        network = self._net
        geoms = self._net.geometries()
        temp = sp.zeros((self._net.Np,))
        for item in geoms:
            ind = network['pore.'+item]
            temp[ind] = temp[ind] + 1
        health = {}
        health['overlaps'] = sp.where(temp>1)[0].tolist()
        health['undefined'] = sp.where(temp==0)[0].tolist()
        return health

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)

