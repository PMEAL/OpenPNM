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
    To invoke the actual generation it is necessary to run the `generate` method.

    Parameters
    ----------
    name : string
        A unique name for the network

    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)
        
    loggername : string
        Overwrite the name of the logger, which defaults to the class name

    Examples
    --------
    >>> pn = OpenPNM.Network.Cubic(name='test_cubic').generate(lattice_spacing=[1],divisions=[5,5,5])

    """

    def __init__(self, **kwargs):
        super(Cubic, self).__init__(**kwargs)
        self._logger.debug(self.__class__.__name__+": Execute constructor")

    def generate(self,**params):
        '''
        Invokes the creation of a Cubic network

        Parameters
        ----------
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

        '''
        self._logger.info(sys._getframe().f_code.co_name+": Start of network topology generation")
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
        self._add_labels()
        self._add_boundaries()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")
        return self

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
        pore_coords = Lc/2+Lc*np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=np.float).T
        self.set_pore_data(prop='coords',data=pore_coords)
        self.set_pore_info(label='all',locations=np.ones_like(ind))
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
        self.set_throat_data(prop='connections',data=connections)      
        self.set_throat_info(label='all',locations=np.ones_like(tpore1))
        self._logger.debug(sys._getframe().f_code.co_name+": End of throat creation")
        
    def _add_labels(self):
        r'''
        Documentation for this method is being updated, we are sorry for the inconvenience.
        '''
        self._logger.info(sys._getframe().f_code.co_name+": Applying labels")
        coords = self.get_pore_data(prop='coords')
        self.set_pore_info(label='front',locations=coords[:,0]<=self._Lc)
        self.set_pore_info(label='left',locations=coords[:,1]<=self._Lc)
        self.set_pore_info(label='bottom',locations=coords[:,2]<=self._Lc)
        self.set_pore_info(label='back',locations=coords[:,0]>=(self._Lc*(self._Nx-1)))
        self.set_pore_info(label='right',locations=coords[:,1]>=(self._Lc*(self._Ny-1)))
        self.set_pore_info(label='top',locations=coords[:,2]>=(self._Lc*(self._Nz-1)))
        self.set_pore_info(label='internal',locations=self.get_pore_indices())
        #Add throat labels based IF both throat's neighbors have label in common
        for item in ['top','bottom','left','right','front','back']:
            ps = self.get_pore_indices(item)
            ts = self.find_neighbor_throats(ps)
            ps = self.find_connected_pores(ts)
            ps0 = self.get_pore_info(label=item)[ps[:,0]]
            ps1 = self.get_pore_info(label=item)[ps[:,1]]
            ts = ts[ps1*ps0]
            self.set_throat_info(label=item,locations=ts)
        self.set_throat_info(label='internal',locations=self.get_throat_indices())
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
            ps = self.get_pore_indices(labels=[label,'internal'],mode='intersection')
            self.clone_pores(pnums=ps,apply_label=[label,'boundary']) 
            #Translate cloned pores
            ind = self.get_pore_indices(labels=[label,'boundary'],mode='intersection')
            coords = self.get_pore_data(prop='coords',locations=ind) 
            coords = coords*scale[label] + offset[label]
            self.set_pore_data(prop='coords', locations=ind, data=coords)

if __name__ == '__main__':
    pn = OpenPNM.Network.Cubic(name='cubic_1',loglevel=10).generate(lattice_spacing=[1.0],domain_size=[3,3,3])
    print(pn.name)
