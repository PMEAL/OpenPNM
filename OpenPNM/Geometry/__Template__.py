"""
module __Template__: Generate cubic networks from domain templates
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry/__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.stats as spst
import scipy.ndimage as spim
from time import clock

from __GenericGeometry__ import GenericGeometry

class Template(GenericGeometry):
    r"""
    Template - Class to create a cubic network with an arbitrary domain shape defined by a supplied template
    
    Parameters
    ----------
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
        
    Examples
    --------
    >>> print 'none yet'
    
    """
    
    def __init__(self, **kwargs):

        super(Template,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        self._logger.info("Import template containing custom network shape")
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")        

        #Instantiate object
        self._net=OpenPNM.Network.GenericNetwork()   
    
    def _generate_setup(self, **params):
        r"""
        
        """
        self._template = params['template']
        self._Lc = params['lattice_spacing']
        if np.ndim(self._template)==3:
            [self._Nx, self._Ny, self._Nz] = np.shape(self._template)
        else:
            [self._Nx, self._Ny] = np.shape(self._template)
            self._Nz = 1
        
    def _generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        
        self._logger.info("generate_pores: Create specified number of pores")
        Lc = self._Lc
        
        #Find non-zero elements in image
        template = self._template
        Np = np.sum(template)
        img_ind = np.ravel_multi_index(np.nonzero(template), dims=np.shape(template), order='F')
        self._net.pore_properties['voxel_index'] = img_ind
        
        #This voxel_to_pore map is messy but works
        temp = np.prod(np.shape(template))*np.ones(np.prod(np.shape(template),),dtype=np.int32)
        temp[img_ind] = np.r_[0:np.size(img_ind)]
        self._voxel_to_pore_map = temp
        
        self._net.pore_properties['coords'] = Lc*(0.5 + np.transpose(np.nonzero(template)))
        self._net.pore_properties['type']= np.zeros((Np,),dtype=np.int8)
        self._net.pore_properties['numbering'] = np.arange(0,Np,dtype=np.int32)
        
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
        img_ind = self._net.pore_properties['voxel_index']
        temp0 = sp.in1d(connections[:,0],img_ind)
        temp1 = sp.in1d(connections[:,1],img_ind)
        tind = temp0*temp1
        connections = connections[tind]
               
        #Need a cleaner way to do this other than voxel_to_pore map...figure out later
        self._net.throat_properties['connections'] = self._voxel_to_pore_map[connections]
        self._net.throat_properties['type'] = np.zeros(np.sum(tind))
        self._net.throat_properties['numbering'] = np.arange(0,np.sum(tind))
        self._logger.debug("generate_throats: End of method")
        
    def _add_boundares(self):
        r"""
        TO DO: Impliment some sort of boundary pore finding
        """
        self._logger.debug("add_boundaries: Nothing yet")
        
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
            net.pore_properties[prop_name] = np.zeros(net.get_num_pores(),dtype=template.dtype)
        
        pore_prop = sp.ravel(template)[net.pore_properties['voxel_index']]
        net.pore_properties[prop_name] = pore_prop
        
        self._logger.debug("add_pore_prop_from_template: End of method")
        
        return pore_prop
        
if __name__ == '__main__':
    test=Template(loggername='TestTemplate')
    test.generate()
    
    
    
    
    
    
    
    