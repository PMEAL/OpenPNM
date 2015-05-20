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
import matplotlib.pyplot as plt
from scipy.io import savemat
from scipy import ndimage
from scipy.stats import itemfreq
logger = logging.getLogger(__name__)


class Voronoi(GenericGeometry):
    r"""
    Voronoi subclass of GenericGeometry.

    Parameters
    ----------

    """

    def __init__(self, fibre_rad=3e-06, load_gen='gen', **kwargs):
        super().__init__(**kwargs)
        self._fibre_rad = fibre_rad
        if load_gen == 'gen':
            self._generate()

    def _generate(self):
        # Set all the required models
        self.models.add(propname='pore.vertices',
                        model=gm.pore_vertices.voronoi)
        self.models.add(propname='throat.vertices',
                        model=gm.throat_vertices.voronoi)
        self.models.add(propname='pore.volume',
                        model=gm.pore_volume.in_hull_volume,
                        fibre_rad=self._fibre_rad)
        self.models.add(propname='throat.normal',
                        model=gm.throat_normal.voronoi)
        self.models.add(propname='throat.offset_vertices',
                        model=gm.throat_offset_vertices.distance_transform,
                        offset=self._fibre_rad,
                        set_dependent=True)
        self.models.add(propname='pore.seed',
                        model=gm.pore_misc.random,
                        seed=self._seed)
        self.models.add(propname='throat.seed',
                        model=gm.throat_misc.neighbor,
                        pore_prop='pore.seed',
                        mode='min')
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
                        const=self._fibre_rad*2)
        self.models.add(propname='throat.volume',
                        model=gm.throat_volume.extrusion)
        self.models.add(propname='throat.surface_area',
                        model=gm.throat_surface_area.extrusion)
        self.models.add(propname='throat.c2c',
                        model=gm.throat_length.voronoi)

    def make_fibre_image(self, fibre_rad=3.5e-6, vox_len=1e-6, add_boundary=False):
        r"""
        If the voronoi voxel method was implemented to calculate pore volumes
        an image of the fibre space has already been calculated and stored on
        the geometry. If not generate it
        """

        if (self._fibre_image is None) or \
           (add_boundary and self.fibre_image_boundary is None):
            self._fibre_image = gm.pore_volume._get_fibre_image(self._net,
                                                                self.pores(),
                                                                vox_len,
                                                                fibre_rad,
                                                                add_boundary)

    def export_fibre_image(self, mat_file='OpenPNMFibres'):
        r"""
        If the voronoi voxel method was implemented to calculate pore volumes
        an image of the fibre space has already been calculated and stored on
        the geometry. If not generate it
        """
        if self._fibre_image is None:
            logger.warning('Fibre image must be generated first')
            return
        matlab_dict = {"fibres": self._fibre_image}
        savemat(mat_file, matlab_dict, format='5', long_field_names=True)

    def get_fibre_slice(self, plane=None, index=None):
        r"""
        Plot an image of a slice through the fibre image
        plane contains percentage values of the length of the image in each axis
        """
        if plane is not None and index is not None:
            logger.warning('Please provide either a plane array or index array')
            return
        if self._fibre_image is None:
            self.make_fibre_image()

        if plane is not None:
            if 'array' not in plane.__class__.__name__:
                plane = sp.asarray(plane)
            if sp.sum(plane == 0) != 2:
                logger.warning('Plane argument must have two zero valued \
                                elements to produce a planar slice')
                return
            l = sp.asarray(sp.shape(self._fibre_image))
            s = sp.around(plane*l).astype(int)
        elif index is not None:
            if 'array' not in index.__class__.__name__:
                index = sp.asarray(index)
            if sp.sum(index == 0) != 2:
                logger.warning('Index argument must have two zero valued \
                                elements to produce a planar slice')
                return
            if 'int' not in str(index.dtype):
                index = sp.around(index).astype(int)
            s = index

        if s[0] != 0:
            slice_image = self._fibre_image[s[0], :, :]
        elif s[1] != 0:
            slice_image = self._fibre_image[:, s[1], :]
        else:
            slice_image = self._fibre_image[:, :, s[2]]

        return slice_image

    def plot_fibre_slice(self, plane=None, index=None):
        r'''
        Plot one slice from the fibre image
        '''
        slice_image = self.get_fibre_slice(plane, index)
        plt.figure()
        plt.imshow(slice_image, cmap='Greys', interpolation='nearest')

    def plot_porosity_profile(self):
        r'''
        Return a porosity profile in all orthogonal directions by summing
        the voxel volumes in consectutive slices.
        '''
        if self._fibre_image is None:
            self.make_fibre_image()
        l = sp.asarray(sp.shape(self._fibre_image))
        px = sp.zeros(l[0])
        py = sp.zeros(l[1])
        pz = sp.zeros(l[2])

        for x in sp.arange(l[0]):
            px[x] = sp.sum(self._fibre_image[x, :, :]) \
                    / sp.size(self._fibre_image[x, :, :])
        for y in sp.arange(l[1]):
            py[y] = sp.sum(self._fibre_image[:, y, :]) \
                    / sp.size(self._fibre_image[:, y, :])
        for z in sp.arange(l[2]):
            pz[z] = sp.sum(self._fibre_image[:, :, z]) \
                    / sp.size(self._fibre_image[:, :, z])

        fig = plt.figure()
        ax = fig.gca()
        plots = []
        plots.append(plt.plot(sp.arange(l[0])/l[0], px, label='x'))
        plots.append(plt.plot(sp.arange(l[1])/l[1], py, label='y'))
        plots.append(plt.plot(sp.arange(l[2])/l[2], pz, label='z'))
        plt.xlabel('Normalized Distance')
        plt.ylabel('Porosity')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=1)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    def _random_array_split(self, array, n):
        '''
        Use array split but shuffle chunks so that bigger and smaller ones are
        randomly distributed.
        '''
        chunks = np.array_split(array, n)
        lengths = np.asarray([len(chunk) for chunk in chunks])
        np.random.shuffle(lengths)
        start = 0
        output = []
        for l in lengths:
            output.append(array[start:start+l])
            start += l
        return output

    def compress_slice(self, slice_image, compression=0.1, plot=False):
        r'''
        compress_slice by shifting fibre pixels
        compression represents the percentage reduction in height
        '''
        shape = sp.shape(slice_image)
        # number of pixels to remove from pore space
        n = sp.ceil(compression*shape[1]).astype(int)
        # initiate everything as pores
        compressed_image = sp.ones([shape[0], shape[1]-n])
        # cycle through columns in the image
        for i in sp.arange(shape[0]):
            column = slice_image[i, :].copy()
            fi = sp.arange(len(column))[column == 0]
            cfi = np.around((1-compression)*fi, 0).astype(int)
            # don't allow fibres to overlap
            freq = itemfreq(cfi)
            # even things out a bit by shifting up and down alternatively
            up_down = 1
            for index, f in freq:
                if f > 1:
                    shift = 1
                    while shift < 20 and f > 1:
                        # check whether next index available
                        test_index = index+(up_down*shift)
                        if (test_index not in cfi) and \
                           (test_index > -1) and \
                           (test_index <= shape[1]-n):
                            cfi = np.hstack((cfi, test_index))
                            f -= 1
                        else:
                            shift += 1
                    up_down *= -1

            column.fill(1)
            column[cfi] = 0
            compressed_image[i, :] = column[:shape[1]-n]

        if plot:
            plt.figure()
            plt.imshow(slice_image, cmap='Greys', interpolation='nearest')
            plt.figure()
            plt.imshow(compressed_image, cmap='Greys', interpolation='nearest')
            print("Fibre voxels before compression: " +
                  str(sp.sum(slice_image == 0)))
            print("Fibre voxels after compression: " +
                  str(sp.sum(compressed_image == 0)))

        return compressed_image

    def compress_fibre_image(self, compression=None):
        r'''
        Compress the fibre image slice by slice
        '''
        fibre_image_shape = np.shape(self._fibre_image)
        num_slice = fibre_image_shape[0]
        compressed_image = np.ones_like(self._fibre_image, dtype=np.uint8)
        for slice_id in sp.arange(1, num_slice):
            logger.info("Compressing Slice "+str(slice_id+1)+" of "+str(num_slice))
            cs = self.compress_slice(slice_image=self._fibre_image[slice_id, :, :],
                                     compression=compression)
            cshape = np.shape(cs)
            compressed_image[slice_id, :cshape[0], :cshape[1]] = cs

        self._compressed_fibre_image = compressed_image[:, :, :cshape[1]]

    def compress_geometry(self, factor=None, preserve_fibres=False):
        r'''
        Adjust the vertices and recalculate geometry. Save fibre voxels before
        and after then put them back into the image to preserve fibre volume.
        Shape will not be conserved. Also make adjustments to the pore and throat
        properties given an approximate volume change from adding the missing fibre
        voxels back in
        '''
        if factor is None:
            logger.warning('Please supply a compression factor in the form \
                            [1,1,CR], with CR < 1')
            return
        if sp.any(sp.asarray(factor) > 1):
            logger.warning('The supplied compression factor is greater than 1, \
                            the method is not tested for stretching')
            return
        # uncompressed number of fibre voxels in each pore
        fvu = self["pore.fibre_voxels"]
        r1 = self["pore.diameter"]/2
        # Most likely boundary pores - prevents divide by zero (vol change zero)
        r1[r1 == 0] = 1
        vo.scale(network=self._net, scale_factor=factor, preserve_vol=False)
        self.models.regenerate()

        if preserve_fibres:
            # compressed number of fibre voxels in each pore
            fvc = self["pore.fibre_voxels"]
            # amount to adjust pore volumes by
            # (based on 1 micron cubed voxels for volume calc)
            vol_diff = (fvu-fvc)*1e-18
            # don't account for positive volume changes
            vol_diff[vol_diff < 0] = 0
            pv1 = self["pore.volume"].copy()
            self["pore.volume"] -= vol_diff
            self["pore.volume"][self["pore.volume"] < 0.0] = 0.0
            pv2 = self["pore.volume"].copy()
            "Now need to adjust the pore diameters"
            from scipy.special import cbrt
            rdiff = cbrt(3*np.abs(vol_diff)/(4*sp.pi))

            self["pore.diameter"] -= 2*rdiff*sp.sign(vol_diff)
            "Now as a crude approximation adjust all the throat areas and diameters"
            "by the same ratio as the increase in a spherical pore surface area"
            spd = np.ones(len(fvu))
            spd[fvu > 0] = (pv2[fvu > 0]/pv1[fvu > 0])**(2/3)
            spd[spd > 1.0] = 1.0
            tconns = self._net["throat.conns"][self.map_throats(self._net,
                                                                self.throats())]
            # Need to work out the average volume change for the two pores
            # connected by each throat Boundary pores will be connected to a
            # throat outside this geometry if there are multiple geoms so get mapping
            mapping = self._net.map_pores(self, self._net.pores(),
                                          return_mapping=True)
            source = list(mapping['source'])
            target = list(mapping['target'])
            ta_diff_avg = np.ones(len(tconns))
            for i in np.arange(len(tconns)):
                np1, np2 = tconns[i]
                if np1 in source and np2 in source:
                    gp1 = target[source.index(np1)]
                    gp2 = target[source.index(np2)]
                    ta_diff_avg[i] = (spd[gp1] + spd[gp2]) / 2
                elif np1 in source:
                    gp1 = target[source.index(np1)]
                    ta_diff_avg[i] = spd[gp1]
                elif np2 in source:
                    gp2 = target[source.index(np2)]
                    ta_diff_avg[i] = spd[gp2]
            self["throat.area"] *= ta_diff_avg
            self["throat.area"][self["throat.area"] < 0] = 0
            self["throat.diameter"] = 2*sp.sqrt(self["throat.area"]/sp.pi)
            self["throat.indiameter"] *= sp.sqrt(ta_diff_avg)
        # Remove pores with zero throats
        self._net.trim_occluded_throats()

if __name__ == '__main__':
    import OpenPNM
    pn = OpenPNM.Network.Delaunay(name='test_net')
    pn.generate(num_pores=100, domain_size=[0.0001, 0.0001, 0.0001])
    pn.add_boundaries()
    test = OpenPNM.Geometry.Voronoi(pores=pn.Ps, throats=pn.Ts, network=pn)
    pn.regenerate_geometries()
