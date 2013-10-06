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
    Template - Class to create a basic cubic network
    
    Parameters
    ----------
    
    domain_size : list with 3 float elements 
        Shape of the cube [Lx,Ly,Lz]
    lattice_spacing : list of three floats
        Spacing between pore centers in each spatial directions
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
        
    Examples
    --------
    
    .. plot::
        
       import pylab as pl
       import OpenPNM
       gen = OpenPNM.Geometry.Cubic()
       net = gen.generate()
       pl.spy(net._adjmatrix)
       pl.show()
    
    TODO:
        - Check for 3D shape
        - Check for correct divisions
    """
    
    def __init__(self,  image_domain = [],
                        image_diameter = [],
                        lattice_spacing = [],
                        **kwargs):

        super(Custom,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        self._logger.info("Import image containing custom network shape")
        
        self._net_img = image_domain
        self._Lc = lattice_spacing
        if np.ndim(image_domain)==3:
            [self._Nx, self._Ny, self._Nz] = np.shape(image_domain)
        self._template = params['template']
        if np.ndim(params['template'])==3:
            [self._Nx, self._Ny, self._Nz] = np.shape(params['template'])
        else:
            [self._Nx, self._Ny] = np.shape(params['template'])
            self._Nz = 1
        Np = self._Nx*self._Ny*self._Nz
        Nt = 3*Np - self._Nx*self._Ny - self._Nx*self._Nz - self._Ny*self._Nz
        
        #Instantiate object
        self._net=OpenPNM.Network.GenericNetwork(num_pores=Np, num_throats=Nt)
    
    def generate_pores(self):
        r"""
        Generate the pores (coordinates, numbering and types)
        """
        self._logger.info("generate_pores: Create specified number of pores")
        Lc = self._Lc
        
        #Find non-zero elements in image
        img = self._template
        Np = np.sum(img)
        img_ind = np.ravel_multi_index(np.nonzero(img), dims=np.shape(img), order='F')
        
        temp = np.prod(np.shape(img))*np.ones(np.prod(np.shape(img),),dtype=np.int32)
        temp[img_ind] = np.r_[0:np.size(img_ind)]
        self._voxel_to_pore_map = temp
        
        self._net.pore_properties['coords'] = Lc*(0.5 + np.transpose(np.nonzero(img)))
        self._net.pore_properties['type']= np.zeros((Np,),dtype=np.int8)
        self._net.pore_properties['numbering'] = np.arange(0,Np,dtype=np.int32)
        
        self._logger.debug("generate_pores: End of method")
        
    def generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        self._logger.info("generate_throats: Define connections between pores")
        
        img = self._template
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
        tpore1 = np.hstack((tpore1_1,tpore1_2,tpore1_3))
        tpore2 = np.hstack((tpore2_1,tpore2_2,tpore2_3))
        connections = np.vstack((tpore1,tpore2)).T
        connections = connections[np.lexsort((connections[:, 1], connections[:, 0]))]
        
        #Remove throats to non-active pores
        img_ind = np.ravel_multi_index(np.nonzero(img), dims=np.shape(img), order='F')
        temp0 = np.in1d(connections[:,0],img_ind)
        temp1 = np.in1d(connections[:,1],img_ind)
        tind = temp0*temp1
        connections = connections[tind]
        
        self._net.throat_properties['connections'] = self._voxel_to_pore_map[connections]
        self._net.throat_properties['type'] = np.zeros(np.sum(tind))
        self._net.throat_properties['numbering'] = np.arange(0,np.sum(tind))
        self._logger.debug("generate_throats: End of method")
#        
#    def generate_pore_seed(self,img1):
#        generate_pore_prop(im_seed,'seed')
#        
#    def generate_pore_diameter(self,img2):
#        generate_pore_prop(img,'diameter')
    def add_boundares(self):
        r"""
        TO DO: Impliment some sort of boundary pore finding
        """
        self._logger.debug("add_boundaries: Nothing yet")
        
    def generate_pore_property_from_image(self,img,prop_name):
        r"""
        Add pore properties based on image location and value
        """
        self._logger.info("add_pore_prop_from_img: Add pore properties")
        
        
        if prop_name not in self._net.pore_properties.keys():
            self._net.pore_properties[prop_name] = np.zeros(self._net.get_num_pores(),dtype=img.dtype)
            
        img_coord = np.nonzero(img) #Find image locations to add
        img_ind = np.ravel_multi_index(img_coord, dims=np.shape(img), order='F')
        vox_map = self._voxel_to_pore_map[img_ind] #Find pore number mapping
        self._net.pore_properties[prop_name][vox_map] = np.array(img[img_coord],dtype=img.dtype)
        
        self._logger.debug("add_pore_prop_from_img: End of method")
        
        
if __name__ == '__main__':
    test=Template(loggername='TestTemplate')
    test.generate()
    
    
    
    
    
    
    
    