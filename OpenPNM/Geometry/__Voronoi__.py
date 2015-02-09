# -*- coding: utf-8 -*-
"""
===============================================================================
Voronoi --Subclass of GenericGeometry for a standard Geometry created from a
Voronoi Diagram Used with Delaunay Network but could work for others (not tested)
===============================================================================


"""

import scipy as sp
from OpenPNM.Geometry import models as gm
from OpenPNM.Geometry import GenericGeometry
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
import matplotlib.pyplot as plt

class Voronoi(GenericGeometry):
    r"""
    Voronoi subclass of GenericGeometry.

    Parameters
    ----------


    """

    def __init__(self, fibre_rad = 3.5e-06, load_gen='gen', **kwargs):
        r"""
        Initialize
        """
        if int(sp.__version__.split('.')[1]) < 13:
            raise Exception('The installed version of Scipy is too old, Voronoi cannot run')
        super(Voronoi,self).__init__(**kwargs)
        if load_gen == 'gen':
            self._generate(fibre_rad)
        

    def _generate(self,fibre_rad):
        r'''
        Set all the required models
        '''
        self.models.add(propname='pore.vertices',
                        model=gm.pore_vertices.voronoi)
        self.models.add(propname='throat.vertices',
                        model=gm.throat_vertices.voronoi)
        self.models.add(propname='throat.normal',
                        model=gm.throat_normal.voronoi)
        self.models.add(propname='throat.offset_vertices',
                        model=gm.throat_offset_vertices.distance_transform,
                        offset=fibre_rad,
                        set_dependent = True)
        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random,
                        seed=self._seed)
        self.models.add(propname='throat.seed',
                        model=gm.throat_misc.neighbor,
                        pore_prop='pore.seed',
                        mode='min')
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.voronoi)
        self.models.add(propname='pore.diameter',
                        model=gm.pore_diameter.voronoi)
        self.models.add(propname='pore.area',
                        model=gm.pore_area.spherical)          
        self.models.add(propname='throat.diameter',
                        model=gm.throat_diameter.voronoi)
        self.models.add(propname='throat.length',
                        model=gm.throat_length.constant,
                        const=fibre_rad*2)
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.extrusion)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.extrusion)
        self.models.add(propname='throat.c2c',
                       model=gm.throat_length.voronoi)

    def make_fibre_image(self,fibre_rad=3.5e-6,vox_len=1e-6):
        r'''
        If the voronoi voxel method was implemented to calculate pore volumes an image of the fibre space
        has already been calculated and stored on the geometry. If not generate it
        '''
        import OpenPNM.Geometry.models as gm
        if self._fibre_image is None:
            self._fibre_image = gm.pore_volume._get_fibre_image(self._net,self.pores(),vox_len,fibre_rad)
        
    def get_fibre_slice(self,plane=None,index=None):
        r'''
        Plot an image of a slice through the fibre image
        plane contains percentage values of the length of the image in each axis
        '''
        if plane is not None and index is not None:
            logger.warning("Please provide either a plane array or index array")
            return
        if self._fibre_image is None:
            self.make_fibre_image()
            
        if plane is not None:
            if 'array' not in plane.__class__.__name__:
                plane = sp.asarray(plane)
            if sp.sum(plane==0) != 2:
                logger.warning("Plane argument must have two zero valued elements to produce a planar slice")
                return
            l = sp.asarray(sp.shape(self._fibre_image))
            s = sp.around(plane*l).astype(int)
        elif index is not None:
            if 'array' not in index.__class__.__name__:
                index = sp.asarray(index)
            if sp.sum(index==0) != 2:
                logger.warning("index argument must have two zero valued elements to produce a planar slice")
                return
            if 'int' not in str(index.dtype):
                index = sp.around(index).astype(int)
            s = index
            
        if s[0] != 0:
            slice_image = self._fibre_image[s[0],:,:]
        elif s[1] != 0:
            slice_image = self._fibre_image[:,s[1],:]
        else:
            slice_image = self._fibre_image[:,:,s[2]]
        
        return slice_image
    
    def plot_fibre_slice(self,plane=None,index=None):
        r'''
        Plot one slice from the fibre image
        '''
        slice_image = self.get_fibre_slice(plane,index)
        plt.figure()
        plt.imshow(slice_image,cmap='Greys',  interpolation='nearest')
        
    def plot_porosity_profile(self):
        r'''
        Return a porosity profile in all orthogonal directions by summing the voxel volumes in consectutive slices
        '''
        if self._fibre_image is None:
            self.make_fibre_image()
        l = sp.asarray(sp.shape(self._fibre_image))
        px = sp.zeros(l[0])
        py = sp.zeros(l[1])
        pz = sp.zeros(l[2])
        
        for x in sp.arange(l[0]):
            px[x] = sp.sum(self._fibre_image[x,:,:])/sp.size(self._fibre_image[x,:,:])
        for y in sp.arange(l[1]):
            py[y] = sp.sum(self._fibre_image[:,y,:])/sp.size(self._fibre_image[:,y,:])
        for z in sp.arange(l[2]):
            pz[z] = sp.sum(self._fibre_image[:,:,z])/sp.size(self._fibre_image[:,:,z])
        
        fig = plt.figure()
        ax = fig.gca()
        plots=[]
        plots.append(plt.plot(sp.arange(l[0])/l[0],px,label='x'))
        plots.append(plt.plot(sp.arange(l[1])/l[1],py,label='y'))
        plots.append(plt.plot(sp.arange(l[2])/l[2],pz,label='z'))
        plt.xlabel('Normalized Distance')
        plt.ylabel('Porosity')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels,loc=1)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    def compress_slice_mark1(self,plane=None,index=None,n=10):
        r'''
        compress_slice by removing pores pixels
        '''
        #n = sp.around(1/compression).astype(int) # number of pixels to remove
        slice_image = self.get_fibre_slice(plane,index)
        shape = sp.shape(slice_image)
        compressed_image = sp.zeros([shape[0],shape[1]-n])
        for i in sp.arange(shape[1]):
            column = slice_image[i,:]
            mini_columns = sp.array_split(sp.arange(len(column)),n)# roughly equally sized arrays of column indices
            ps = sp.array([],dtype=int)
            for mini in mini_columns:
                pi = mini[column[mini]==1]# non-zero indices (i.e. pore space)
                fi = mini[column[mini]==0]# zero indices (i.e. fibre space)
                if len(pi > 0):                           
                    #add back in pore pixels missing one from every mini-column           
                    ps = sp.concatenate((ps,fi,pi[:len(pi)-1]))
                else:
                    #destroying fibre - bad but work out later!!!
                    ps = sp.concatenate((ps,fi[:len(fi)-1])) 
            ps.sort()
            compressed_image[i,:]=column[ps]
        plt.figure()
        plt.imshow(slice_image,cmap='Greys',  interpolation='nearest')
        plt.figure()
        plt.imshow(compressed_image,cmap='Greys',  interpolation='nearest')
        print("Fibre voxels before compression: "+str(sp.sum(slice_image==0)))
        print("Fibre voxels after compression: "+str(sp.sum(compressed_image==0)))
        
    def compress_slice(self,plane=None,index=None,compression=0.1):
        r'''
        compress_slice by removing pores pixels
        '''
        #n = sp.around(1/compression).astype(int) # number of pixels to remove
        slice_image = self.get_fibre_slice(plane,index)
        shape = sp.shape(slice_image)
        n = sp.ceil(compression*shape[1]).astype(int) #number of pixels to remove from pore space
        compressed_image = sp.zeros([shape[0],shape[1]-n])#initiate everything as fibre
        #cycle through columns in the image
        for i in sp.arange(shape[1]):
            column = slice_image[i,:]
            pi = sp.arange(len(column))[column==1]# non-zero indices (i.e. pore space)
            fi = sp.arange(len(column))[column==0]# zero indices (i.e. fibre space)
            sp.random.shuffle(pi) # randomize the pore space indices
            if len(pi) > n:
                ps = sp.array([],dtype=int)
                chunks = sp.array_split(pi,n)
                for chunk in chunks:
                    ps = sp.concatenate((ps,chunk[0:len(chunk)-1]))
                ps = sp.concatenate((fi,ps))
                ps.sort()
                compressed_image[i,:len(ps)]=column[ps]
            else:
                #everything fibre
                pass
            
        plt.figure()
        plt.imshow(slice_image,cmap='Greys',  interpolation='nearest')
        plt.figure()
        plt.imshow(compressed_image,cmap='Greys',  interpolation='nearest')
        print("Fibre voxels before compression: "+str(sp.sum(slice_image==0)))
        print("Fibre voxels after compression: "+str(sp.sum(compressed_image==0)))
        
if __name__ == '__main__':
    import OpenPNM
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001])
    pn.add_boundaries()
    test = OpenPNM.Geometry.Voronoi(pores=pn.Ps,throats=pn.Ts,network=pn)
    pn.regenerate_geometries()