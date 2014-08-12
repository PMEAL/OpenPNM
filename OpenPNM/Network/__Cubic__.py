"""
module __Cubic__: Generate simple cubic networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Topology.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
import numpy as np
import scipy.stats as spst
import scipy.spatial as sptl
import itertools as itr
from OpenPNM.Network import GenericNetwork


class Cubic(GenericNetwork):
    r"""
    This class contains the methods for creating a *Cubic* network topology.  
    The actual generation is carried out during the initialization based upon the 
    passed parameters domain_size, lattice_spacing and divisions.
    
    To invoke the actual generation it is necessary to run the `generate` method.

    Parameters
    ----------
    name : string
        A unique name for the network

    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)
        
    loggername : string
        Overwrite the name of the logger, which defaults to the class name

    domain_size : [float,float,float]
        domain_size = [3.0,3.0,3.0] (default)\n
        Bounding cube for internal pore positions\n
    lattice_spacing : [float]
        lattice_spacing = [1.0] (default)\n
        Distance between pore centers\n
    divisions : [int,int,int]
        divisions = [3,3,3]\n
        Number of internal pores in each dimension.\n
        (Optional input. Replaces one of the above.)\n


    Examples
    --------
    >>> pn = OpenPNM.Network.Cubic()
    >>> pn.generate(lattice_spacing=[1],divisions=[5,5,5],add_boundaries=False)

    """

    def __init__(self, **params):
        super(Cubic, self).__init__(**params)
        self._logger.debug(self.__class__.__name__+": Execute constructor")
        
        self._logger.info(sys._getframe().f_code.co_name+": Start of network topology generation")
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
        self._add_labels()
        if params['add_boundaries']:
            self._add_boundaries()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")
        
    def _generate_setup(self,   domain_size = [],
                                divisions = [],
                                lattice_spacing = [],
                                btype = [0,0,0],
                                **params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")
        #Parse the given network size variables
        self._btype = btype
        if domain_size and lattice_spacing and not divisions:
            self._Lc = np.float(lattice_spacing[0])
            self._Lx = np.float(domain_size[0])
            self._Ly = np.float(domain_size[1])
            self._Lz = np.float(domain_size[2])
            self._Nx = np.int(self._Lx/self._Lc)
            self._Ny = np.int(self._Ly/self._Lc)
            self._Nz = np.int(self._Lz/self._Lc)
        elif divisions and lattice_spacing and not domain_size:
            self._Lc = np.float(lattice_spacing[0])
            self._Nx = np.int(divisions[0])
            self._Ny = np.int(divisions[1])
            self._Nz = np.int(divisions[2])
            self._Lx = np.float(self._Nx*self._Lc)
            self._Ly = np.float(self._Ny*self._Lc)
            self._Lz = np.float(self._Nz*self._Lc)
        elif domain_size and divisions and not lattice_spacing:
            self._Lc = np.min(np.array(domain_size,dtype=np.float)/np.array(divisions,dtype=np.float))
            self._Nx = np.int(divisions[0])
            self._Ny = np.int(divisions[1])
            self._Nz = np.int(divisions[2])
            self._Lx = np.float(self._Nx*self._Lc)
            self._Ly = np.float(self._Ny*self._Lc)
            self._Lz = np.float(self._Nz*self._Lc)
        elif not domain_size and not divisions and not lattice_spacing:
            self._Lc = np.float(1)
            self._Lx = np.float(3)
            self._Ly = np.float(3)
            self._Lz = np.float(3)
            self._Nx = np.int(self._Lx/self._Lc)
            self._Ny = np.int(self._Ly/self._Lc)
            self._Nz = np.int(self._Lz/self._Lc)
        else:
            self._logger.error("Exactly two of domain_size, divisions and lattice_spacing must be given")
            raise Exception('Exactly two of domain_size, divisions and lattice_spacing must be given')

    def _generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        self._logger.info(sys._getframe().f_code.co_name+": Creating specified number of pores")
        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Lc = self._Lc
        Np = Nx*Ny*Nz
        ind = np.arange(0,Np)
        self['pore.coords'] = Lc/2+Lc*np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=np.float).T
        self['pore.all'] = np.ones_like(ind,dtype=bool)
        self._logger.debug(sys._getframe().f_code.co_name+": End of pore creation")

    def _generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        self._logger.info(sys._getframe().f_code.co_name+": Creating specified number of throats")

        Nx = self._Nx
        Ny = self._Ny
        Nz = self._Nz
        Np = Nx*Ny*Nz
        ind = np.arange(0,Np)

        #Generate throats based on pattern of the adjacency matrix
        tpore1_1 = ind[(ind%Nx)<(Nx-1)]
        tpore2_1 = tpore1_1 + 1
        tpore1_2 = ind[(ind%(Nx*Ny))<(Nx*(Ny-1))]
        tpore2_2 = tpore1_2 + Nx
        tpore1_3 = ind[(ind%Np)<(Nx*Ny*(Nz-1))]
        tpore2_3 = tpore1_3 + Nx*Ny
        tpore1 = np.hstack((tpore1_1,tpore1_2,tpore1_3))
        tpore2 = np.hstack((tpore2_1,tpore2_2,tpore2_3))
        connections = np.vstack((tpore1,tpore2)).T
        connections = connections[np.lexsort((connections[:, 1], connections[:, 0]))]
        self['throat.conns'] = connections
        self['throat.all'] = np.ones_like(tpore1,dtype=bool)
        self._logger.debug(sys._getframe().f_code.co_name+": End of throat creation")
        
    def _add_labels(self):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        self._logger.info(sys._getframe().f_code.co_name+": Applying labels")
        coords = self['pore.coords']
        self['pore.front'] = self.tomask(coords[:,0]<=self._Lc)
        self['pore.left'] = self.tomask(coords[:,1]<=self._Lc)
        self['pore.bottom'] = self.tomask(coords[:,2]<=self._Lc)
        self['pore.back'] = self.tomask(coords[:,0]>=(self._Lc*(self._Nx-1)))
        self['pore.right'] = self.tomask(coords[:,1]>=(self._Lc*(self._Ny-1)))
        self['pore.top'] = self.tomask(coords[:,2]>=(self._Lc*(self._Nz-1)))
        for item in ['top','bottom','left','right','front','back']:
            ps = self.pores(item)
            ts = self.find_neighbor_throats(ps)
            ps = self.find_connected_pores(ts)
            ps0 = self['pore.'+item][ps[:,0]]
            ps1 = self['pore.'+item][ps[:,1]]
            ts = ts[ps1*ps0]
            self['throat.'+item] = self.tomask(throats=ts)
        self._logger.debug(sys._getframe().f_code.co_name+": End")
        
    def _add_boundaries(self):
        r'''
        This method uses clone_pore to clone the surface pores (labeled 'left'
        , 'right', etc), then shifts them to the periphery of the domain.
        '''
        
        offset = {}
        offset['front'] = offset['left'] = offset['bottom'] = [0,0,0]
        offset['back']  = [self._Lx,0,0]
        offset['right'] = [0,self._Ly,0]
        offset['top']   = [0,0,self._Lz]
        
        scale = {}
        scale['front']  = scale['back']  = [0,1,1]
        scale['left']   = scale['right'] = [1,0,1]
        scale['bottom'] = scale['top']   = [1,1,0]
        
        for label in ['front','back','left','right','bottom','top']:
            ps = self.pores(label)
            self.clone(pores=ps,apply_label=['boundary',label+'_face',label])
            #Translate cloned pores
            ind = self.pores(label+'_face')
            coords = self['pore.coords'][ind] 
            coords = coords*scale[label] + offset[label]
            self['pore.coords'][ind] = coords
            
    def toarray(self,data=None):
        r'''
        Returns a cubic array with each element representing a pore, containing
        specified data values
        
        Parameters
        ----------
        data : array_like
            An Np-long array containing pore property values to enter in each
            array element
        '''
        if data == None:
            data = sp.ones((self.Np,))
        temp = sp.reshape(data,(self._Nx,self._Ny,self._Nz))
        return temp
        
        
if __name__ == '__main__':
    pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=10).generate(lattice_spacing=[1.0],domain_size=[3,3,3])
    print(pn.name)
