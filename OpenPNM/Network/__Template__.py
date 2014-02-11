"""
module __Template__: Generate cubic networks from domain templates
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry/__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

import scipy as sp
import numpy as np
from OpenPNM.Network.__GenericNetwork__ import GenericNetwork

class Template(GenericNetwork):
    r"""
    This class contains the methods to create a cubic network with an arbitrary 
    domain shape defined by a supplied template.
    
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
    >>> img = sp.ones((30,30,30),dtype=int)
    >>> pn = OpenPNM.Network.Template(name='template_1').generate(template=img,lattice_spacing=0.001)
    >>> pn.num_pores()
    27000
    >>> pn.num_throats()
    78300
    """

    def __init__(self, **kwargs):
        super(Template,self).__init__(**kwargs)
        self._logger.debug(self.__class__.__name__,": ","Execute constructor")
    
    def generate(self, **params):
        r'''
        template : array_like, boolean
            An image containing 1's where each pore should be located, and 0's elsewhere.
            This image can be 1,2 or 3D.
        lattice_spacing : float
            The lattice constant for the network, used to scale distance between pores.
        '''
        self._logger.info(sys._getframe().f_code.co_name+": Start of network topology generation")
        self._generate_setup(**params)
        self._generate_pores()
        self._generate_throats()
        self._add_boundaries()
#        self._add_labels()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")
        return self

    def _generate_setup(self, **params):
        r"""

        """
        self._template = params['template']
        if params['lattice_spacing']:
            self._Lc = params['lattice_spacing']
        else:
            self._logger.error("lattice_spacing not specified")
            raise Exception('lattice_spacing not specified')
        if np.ndim(self._template)==3:
            [self._Nx, self._Ny, self._Nz] = np.shape(self._template)
        elif np.ndim(self._template)==2:
            [self._Nx, self._Ny] = np.shape(self._template)
            self._Nz = 1
        elif np.ndim(self._template)==1:
            self._Nx = np.shape(self._template)
            self._Ny = 1
            self._Nz = 1

    def _generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        self._logger.info("generate_pores: Create specified number of pores")
        Lc = self._Lc

        #Find non-zero elements in image
        template = self._template
        Np = np.sum(template>0)
        img_ind = np.ravel_multi_index(np.nonzero(template), dims=np.shape(template), order='F')
        self.set_pore_data(prop='voxel_index',data=img_ind)

        #This voxel_to_pore map is messy but works
        temp = np.prod(np.shape(template))*np.ones(np.prod(np.shape(template),),dtype=np.int32)
        temp[img_ind] = np.r_[0:np.size(img_ind)]
        self._voxel_to_pore_map = temp

        if self._Nz == 1:
            print(np.shape(Lc*(0.5 + np.transpose(np.nonzero(template)))))
            print(np.shape(np.zeros((Np,1))))
            coords = sp.hstack((Lc*(0.5 + np.transpose(np.nonzero(template))),np.zeros((Np,1))))
            self.set_pore_data(prop='coords',data=coords)
        else:
            coords = Lc*(0.5 + np.transpose(np.nonzero(template)))
            self.set_pore_data(prop='coords',data=coords)
        ind = np.arange(0,Np,dtype=np.int32)
        self.set_pore_data(prop='numbering',data=ind)
        self.set_pore_info(prop='all',locations=sp.ones_like(ind))
        self._logger.debug("generate_pores: End of method")

    def _generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        self._logger.info("generate_throats: Define connections between pores")

        template = self._template
        [Nx, Ny, Nz] = [self._Nx, self._Ny, self._Nz]
        Np = Nx*Ny*Nz
        ind = np.arange(0,Np)

        #Generate throats based on pattern of the adjacency matrix
        tpore1_1 = ind[(ind%Nx)<(Nx-1)]
        tpore2_1 = tpore1_1 + 1
        tpore1_2 = ind[(ind%(Nx*Ny))<(Nx*(Ny-1))]
        tpore2_2 = tpore1_2 + Nx
        tpore1_3 = ind[(ind%Np)<(Nx*Ny*(Nz-1))]
        tpore2_3 = tpore1_3 + Nx*Ny
        tpore1 = sp.hstack((tpore1_1,tpore1_2,tpore1_3))
        tpore2 = sp.hstack((tpore2_1,tpore2_2,tpore2_3))
        connections = sp.vstack((tpore1,tpore2)).T
        connections = connections[sp.lexsort((connections[:, 1], connections[:, 0]))]

        #Remove throats to non-active pores
        img_ind = self.get_pore_data(prop='voxel_index')
        temp0 = sp.in1d(connections[:,0],img_ind)
        temp1 = sp.in1d(connections[:,1],img_ind)
        tind = temp0*temp1
        connections = connections[tind]

        #Need a cleaner way to do this other than voxel_to_pore map...figure out later
        self.set_throat_data(prop='connections',data=self._voxel_to_pore_map[connections])
        self.set_pore_data(prop='numbering', data=np.arange(0,np.sum(tind)))
        self.set_throat_info(prop='all',locations=sp.ones_like(tind))
        self._logger.debug("generate_throats: End of method")

    def _add_boundaries(self):
        r"""
        TO DO: Impliment some sort of boundary pore finding
        """
        self._logger.debug("add_boundaries: Not implemented")


    def add_pore_property_from_template(self,net,template,prop_name):
        r"""
        Add pore properties based on value stored at each location in the template array

        Parameters
        ----------
        net : Pore Network Object
            The network to which the pore property should be added

        template : array_like
            The template array containing the pore property values at the desired locations

        prop_name : string
            The name of the pore property being added

        Returns
        -------
        pore_prop : array_like, optional
            The pore_property array extracted from the template image can be received if desire, but it is also written to the object's pore_property dictionary under the key prop_name.

        Notes
        -----
        This method can lead to troubles if not executed in the right order.  For instance, if throat sizes are assigned during the generation stage based on neighboring pores sizes, then rewriting pore sizes with this method could invalidate the throat sizes.  Ideally, this method should be called by the generate_pore_sizes() step during the generate process to avoid the issue.  Alternatively, an update_throat_sizes() method could be written and called subsequent to calling this method.

        """
        self._logger.info("add_pore_prop_from_template: Add pore properties")
        if prop_name not in net.pore_properties.keys():
            net.pore_properties[prop_name] = np.zeros(net.num_pores(),dtype=template.dtype)
        pore_prop = sp.ravel(template)[net.pore_properties['voxel_index']]
        net.pore_properties[prop_name] = pore_prop
        self._logger.debug("add_pore_prop_from_template: End of method")
        return pore_prop

if __name__ == '__main__':
    tmplt = sp.ones((30,30,30),dtype=int)
    pn = OpenPNM.Network.Template(name='template_1').generate(template=tmplt,lattice_spacing=0.001)
    print(pn.name)







