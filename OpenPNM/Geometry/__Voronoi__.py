# -*- coding: utf-8 -*-
"""
===============================================================================
Voronoi --Subclass of GenericGeometry for a standard Geometry created from a
Voronoi Diagram Used with Delaunay Network but could work for others (not tested)
===============================================================================


"""

import scipy as sp
import numpy as np
from OpenPNM.Geometry import models as gm
import OpenPNM.Utilities.vertexops as vo
from OpenPNM.Geometry import GenericGeometry
from OpenPNM.Base import logging
logger = logging.getLogger(__name__)
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy import ndimage
from scipy.stats import itemfreq

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
                        model=gm.pore_volume.voronoi_vox,fibre_rad=fibre_rad)
        self.models.add(propname='pore.centroid',
                        model=gm.pore_centroid.voronoi)
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

    def make_fibre_image(self,fibre_rad=3.5e-6,vox_len=1e-6,add_boundary=False):
        r'''
        If the voronoi voxel method was implemented to calculate pore volumes an image of the fibre space
        has already been calculated and stored on the geometry. If not generate it
        '''
        
        if (self._fibre_image is None) or (add_boundary and self.fibre_image_boundary is None):
            self._fibre_image = gm.pore_volume._get_fibre_image(self._net,self.pores(),vox_len,fibre_rad,add_boundary)
    
    def export_fibre_image(self,mat_file='OpenPNMFibres'):
        r'''
        If the voronoi voxel method was implemented to calculate pore volumes an image of the fibre space
        has already been calculated and stored on the geometry. If not generate it
        '''
        if self._fibre_image is None:
            logger.warning("Fibre image must be generated first")
            return
        matlab_dict = {"fibres":self._fibre_image}
        savemat(mat_file,matlab_dict,format='5',long_field_names=True)
        
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
        
    def compress_geometry(self,factor=None):
        r'''
        Adjust the vertices and recalculate geometry. Save fibre voxels before and after then put them 
        back into the image to preserve fibre volume. Shape will not be conserved. Also make adjustments 
        to the pore and throat properties given an approximate volume change from adding the missing fibre
        voxels back in
        '''
        if factor is None:
            logger.warning("Please supply a compression factor in the form [1,1,CR], with CR < 1")
            return
        if sp.any(sp.asarray(factor)>1):
            logger.warning("The supplied compression factor is greater than 1, the method is not tested for stretching")
            return
        fvu = self["pore.fibre_voxels"] # uncompressed number of fibre voxels in each pore
        vo.scale(network=self._net,scale_factor=factor,preserve_vol=False)
        self.models.regenerate()
        fvc = self["pore.fibre_voxels"] # compressed number of fibre voxels in each pore
        vox_diff = fvu-fvc # Number of fibre voxels to put back into each pore
        n = sp.sum(vox_diff) # Total number of voxels to put back into image
        vol_diff = (fvu-fvc)*1e-18 # amount to adjust pore volumes by
        self["pore.volume"] -= vol_diff
        "Now need to adjust the pore diameters"
        from scipy.special import cbrt
        rdiff = cbrt(3*np.abs(vol_diff)/(4*sp.pi))
        r1 = self["pore.diameter"]/2
        self["pore.diameter"] -= 2*rdiff*sp.sign(vol_diff)
        "Now as a crude approximation adjust all the throat areas and diameters"
        "by the same ratio as the increase in a spherical pore surface area"
        spd = 2*rdiff/r1 + (rdiff/r1)**2 # surface-area percentage difference
        tconns = self._net["throat.conns"][self.map_throats(self._net,self.throats())]
        "Need to work out the average volume change for the two pores connected by each throat"
        "Boundary pores will be connected to a throat outside this geometry if there are multiple geoms so get mapping"
        mapping = self._net.map_pores(self,self._net.pores(),return_mapping=True)
        source = list(mapping['source'])
        target = list(mapping['target'])
        ta_diff_avg = np.zeros(len(tconns))
        for i in np.arange(len(tconns)):
            np1,np2 = tconns[i]
            if np1 in source and np2 in source:
                gp1 = target[source.index(np1)]
                gp2 = target[source.index(np2)]
                ta_diff_avg[i] = (spd[gp1]*sp.sign(vol_diff[gp1])+spd[gp2]*sp.sign(vol_diff[gp2]))/2
            elif np1 in source:
                gp1 = target[source.index(np1)]
                ta_diff_avg[i] = spd[gp1]*sp.sign(vol_diff[gp1])
            elif np2 in mapping['source']:
                gp2 = target[source.index(np2)]
                ta_diff_avg[i] = spd[gp2]*sp.sign(vol_diff[gp2])
                
        self["throat.area"] *= 1+ta_diff_avg
        self["throat.area"][self["throat.area"]<0]=0
        self["throat.diameter"] = 2*sp.sqrt(self["throat.area"]/sp.pi)
        "Now the recreated fibre image will have the wrong number of fibre volumes so add them back in at semi-random"
        b_img = self._fibre_image_boundary #distance transform image to identify boundary layer
        "Cycle through to find distance corresponding roughly to where new fibre layer should be"
        for l in itemfreq(b_img)[1:,0]:
            sel = (b_img == l)
            b = sp.sum(sel) #number of voxels in current boundary layer
            if b < n:
                self._fibre_image[sel]=0 # expand fibre space to fill whole boundary layer
                n -= b # adjust n to be remainder
            else:
                if n > 0:
                    m = 0
                    for i in sp.nditer(sel,op_flags=['readwrite']): #only keep n True values in the selection
                        if i == True and m >= n:
                            i[...] = False
                            m +=1
                        elif i == True:
                            m +=1
                    self._fibre_image[sel]=0 # expand fibre space to fill first n boundary layer voxels
                    break
                else:
                    break
        if sp.sum(self._fibre_image==0) != sp.sum(fvu):
            print("Something went wrong with compression")
if __name__ == '__main__':
    import OpenPNM
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001,0.0001,0.0001])
    pn.add_boundaries()
    test = OpenPNM.Geometry.Voronoi(pores=pn.Ps,throats=pn.Ts,network=pn)
    pn.regenerate_geometries()