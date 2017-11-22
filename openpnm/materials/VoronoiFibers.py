import scipy as sp
import numpy as np
import openpnm.utils.vertexops as vo
import scipy.sparse as sprs
import scipy.spatial as sptl
import matplotlib.pyplot as plt
from scipy.io import savemat
from openpnm import topotools
from openpnm.network import DelaunayVoronoiDual
from openpnm.core import logging
from openpnm.geometry import models as gm
from openpnm.geometry import GenericGeometry
logger = logging.getLogger(__name__)


class VoronoiFibers(DelaunayVoronoiDual):
    r"""

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        DelaunayGeometry(network=self, pores=self.pores('delaunay'),
                         throats=self.throats('delaunay'),
                         name=self.name+'_del')


class DelaunayGeometry(GenericGeometry):
    r"""
    Voronoi subclass of GenericGeometry.

    Parameters
    ----------
    name : string
        A unique name for the network

    fibre_rad: float
        Fibre radius to apply to Voronoi edges when calculating pore and throat
        sizes

    voxel_vol : boolean
        Determines whether to calculate pore volumes by creating a voxel image
        or to use the offset vertices of the throats. Voxel method is slower
        and may run into memory issues but is more accurate and allows manipulation
        of the image. N.B. many of the class methods are dependent on the voxel
        image.
    """

    def __init__(self, network, fibre_rad=3e-06, voxel_vol=True, **kwargs):
        super().__init__(network=network, **kwargs)
        self._fibre_rad = fibre_rad
        self._voxel_vol = voxel_vol
        if 'vox_len' in kwargs.keys():
            self._vox_len = kwargs['vox_len']
        else:
            self._vox_len = 1e-6

        # Set all the required models
        vertices = network.find_pore_hulls()
        coords = [network['pore.coords'][p] for p in vertices]
        self['pore.vertices'] = coords
        vertices = network.find_throat_facets()
        coords = [network['pore.coords'][t] for t in vertices]
        self['throat.vertices'] = coords

        self.add_model(propname='throat.normal',
                       model=gm.throat_normal.voronoi)
        self.add_model(propname='throat.offset_vertices',
                       model=gm.throat_offset_vertices.distance_transform,
                       offset=self._fibre_rad,
                       set_dependent=True)

        topotools.trim_occluded_throats(network=network, mask=self.name)

#        if self._voxel_vol:
#            self.add_model(propname='pore.volume',
#                           model=gm.pore_volume.in_hull_volume,
#                           fibre_rad=self._fibre_rad,
#                           vox_len=self._vox_len)
#        else:
#            self.add_model(propname='pore.volume',
#                           model=gm.pore_volume.voronoi)
#        self.add_model(propname='throat.shape_factor',
#                       model=gm.throat_shape_factor.compactness)
#        self.add_model(propname='pore.seed',
#                       model=gm.pore_misc.random)
#        self.add_model(propname='throat.seed',
#                       model=gm.throat_misc.neighbor,
#                       pore_prop='pore.seed',
#                       mode='min')
#        self.add_model(propname='pore.centroid',
#                       model=gm.pore_centroid.voronoi)
#        self.add_model(propname='pore.diameter',
#                       model=gm.pore_diameter.equivalent_sphere)
#        self.add_model(propname='pore.indiameter',
#                       model=gm.pore_diameter.centroids)
#        self.add_model(propname='pore.area',
#                       model=gm.pore_area.spherical)
#        self.add_model(propname='throat.diameter',
#                       model=gm.throat_diameter.equivalent_circle)
#        self['throat.volume'] = 0.0
#        self['throat.length'] = self._fibre_rad*2
#        self.add_model(propname='throat.surface_area',
#                       model=gm.throat_surface_area.extrusion)
#        self.add_model(propname='throat.c2c',
#                       model=gm.throat_length.c2c)

    def make_fibre_image(self, fibre_rad=None, vox_len=1e-6):
        r"""
        If the voronoi voxel method was implemented to calculate pore volumes
        an image of the fibre space has already been calculated and stored on
        the geometry. If not generate it

        Parameters
        ----------
        fibre_rad : float
        Fibre radius to apply to Voronoi edges when calculating pore and throat
        sizes

        vox_len : float
        Length of voxel edge when dividing domain
        """

        if hasattr(self, '_fibre_image'):
            logger.info('fibre image already created')
            return
        else:
            if fibre_rad is None:
                fibre_rad = self._fibre_rad
            fibre_rad /= vox_len
            self._fibre_image = gm.pore_volume._get_fibre_image(self._net,
                                                                self.pores(),
                                                                vox_len,
                                                                fibre_rad)

    def _get_fibre_image(self, vox_len, fibre_rad):
        r"""
        Produce image by filling in voxels along throat edges using Bresenham line
        Then performing distance transform on fibre voxels to erode the pore space
        """
        network = self.simulation.network
        cthroats = network.find_neighbor_throats(pores=cpores)

        # Below method copied from geometry model throat.vertices
        # Needed now as network may not have all throats assigned to geometry
        # i.e network['throat.vertices'] could return garbage
        verts = _sp.ndarray(network.num_throats(), dtype=object)
        for i in range(len(verts)):
            verts[i] = _sp.asarray(list(network["throat.vert_index"][i].values()))
        cverts = verts[cthroats]
        [vxmin, vxmax, vymin, vymax, vzmin, vzmax] = _get_vertex_range(cverts)
        # Translate vertices so that minimum occurs at the origin
        for index in range(len(cverts)):
            cverts[index] -= np.array([vxmin, vymin, vzmin])
        # Find new size of image array
        cdomain = np.around(np.array([(vxmax-vxmin),
                                      (vymax-vymin),
                                      (vzmax-vzmin)]), 6)
        logger.info("Creating fibre domain range: " + str(np.around(cdomain, 5)))
        lx = np.int(np.around(cdomain[0]/vox_len)+1)
        ly = np.int(np.around(cdomain[1]/vox_len)+1)
        lz = np.int(np.around(cdomain[2]/vox_len)+1)
        # Try to create all the arrays we will need at total domain size
        try:
            pore_space = np.ones([lx, ly, lz], dtype=np.uint8)
            fibre_space = np.zeros(shape=[lx, ly, lz], dtype=np.uint8)
            dt = np.zeros([lx, ly, lz], dtype=float)
            # Only need one chunk
            cx = cy = cz = 1
            chunk_len = np.max(np.shape(pore_space))
        except:
            logger.info("Domain too large to fit into memory so chunking domain"
                        "to process image, this may take some time")
            # Do chunking
            chunk_len = 100
            if (lx > chunk_len):
                cx = np.ceil(lx/chunk_len).astype(int)
            else:
                cx = 1
            if (ly > chunk_len):
                cy = np.ceil(ly/chunk_len).astype(int)
            else:
                cy = 1
            if (lz > chunk_len):
                cz = np.ceil(lz/chunk_len).astype(int)
            else:
                cz = 1

        # Get image of the fibres
        line_points = bresenham(cverts, vox_len/2)
        line_ints = (np.around((line_points/vox_len), 0)).astype(int)
        for x, y, z in line_ints:
            try:
                pore_space[x][y][z] = 0
            except IndexError:
                logger.warning("Some elements in image processing are out" +
                               "of bounds")

        num_chunks = np.int(cx*cy*cz)
        cnum = 1
        for ci in range(cx):
            for cj in range(cy):
                for ck in range(cz):
                    # Work out chunk range
                    logger.info("Processing Fibre Chunk: "+str(cnum)+" of " +
                                str(num_chunks))
                    cxmin = ci*chunk_len
                    cxmax = np.int(np.ceil((ci+1)*chunk_len + 5*fibre_rad))
                    cymin = cj*chunk_len
                    cymax = np.int(np.ceil((cj+1)*chunk_len + 5*fibre_rad))
                    czmin = ck*chunk_len
                    czmax = np.int(np.ceil((ck+1)*chunk_len + 5*fibre_rad))
                    # Don't overshoot
                    if cxmax > lx:
                        cxmax = lx
                    if cymax > ly:
                        cymax = ly
                    if czmax > lz:
                        czmax = lz
                    dt = ndimage.distance_transform_edt(pore_space[cxmin:cxmax,
                                                                   cymin:cymax,
                                                                   czmin:czmax])
                    fibre_space[cxmin:cxmax,
                                cymin:cymax,
                                czmin:czmax][dt <= fibre_rad] = 0
                    fibre_space[cxmin:cxmax,
                                cymin:cymax,
                                czmin:czmax][dt > fibre_rad] = 1
                    del dt
                    cnum += 1
        del pore_space
        return fibre_space

    def _export_fibre_image(self, mat_file='OpenPNMFibres'):
        r"""
        If the voronoi voxel method was implemented to calculate pore volumes
        an image of the fibre space has already been calculated and stored on
        the geometry.

        Parameters
        ----------
        mat_file : string
        Filename of Matlab file to save fibre image
        """
        if hasattr(self, '_fibre_image') is False:
            logger.warning('This method only works when a fibre image exists, ' +
                           'please run make_fibre_image')
            return
        matlab_dict = {"fibres": self._fibre_image}
        savemat(mat_file, matlab_dict, format='5', long_field_names=True)

    def _get_fibre_slice(self, plane=None, index=None):
        r"""
        Plot an image of a slice through the fibre image
        plane contains percentage values of the length of the image in each axis

        Parameters
        ----------
        plane : array_like
        List of 3 values, [x,y,z], 2 must be zero and the other must be between
        zero and one representing the fraction of the domain to slice along
        the non-zero axis

        index : array_like
        similar to plane but instead of the fraction an index of the image is used
        """
        if hasattr(self, '_fibre_image') is False:
            logger.warning('This method only works when a fibre image exists, ' +
                           'please run make_fibre_image')
            return None
        if plane is None and index is None:
            logger.warning('Please provide either a plane array or index array')
            return None
        if self._fibre_image is None:
            self.make_fibre_image()

        if plane is not None:
            if 'array' not in plane.__class__.__name__:
                plane = sp.asarray(plane)
            if sp.sum(plane == 0) != 2:
                logger.warning('Plane argument must have two zero valued ' +
                               'elements to produce a planar slice')
                return None
            l = sp.asarray(sp.shape(self._fibre_image))
            s = sp.around(plane*l).astype(int)
        elif index is not None:
            if 'array' not in index.__class__.__name__:
                index = sp.asarray(index)
            if sp.sum(index == 0) != 2:
                logger.warning('Index argument must have two zero valued ' +
                               'elements to produce a planar slice')
                return None
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

    def plot_fibre_slice(self, plane=None, index=None, fig=None):
        r"""
        Plot one slice from the fibre image

        Parameters
        ----------
        plane : array_like
        List of 3 values, [x,y,z], 2 must be zero and the other must be between
        zero and one representing the fraction of the domain to slice along
        the non-zero axis

        index : array_like
        similar to plane but instead of the fraction an index of the image is used
        """
        if hasattr(self, '_fibre_image') is False:
            logger.warning('This method only works when a fibre image exists, ' +
                           'please run make_fibre_image')
            return
        slice_image = self._get_fibre_slice(plane, index)
        if slice_image is not None:
            if fig is None:
                plt.figure()
            plt.imshow(slice_image.T, cmap='Greys', origin='lower',
                       interpolation='nearest')

        return fig

    def plot_porosity_profile(self, fig=None):
        r"""
        Return a porosity profile in all orthogonal directions by summing
        the voxel volumes in consectutive slices.
        """
        if hasattr(self, '_fibre_image') is False:
            logger.warning('This method only works when a fibre image exists, ' +
                           'please run make_fibre_image')
            return
        if self._fibre_image is None:
            self.make_fibre_image()
        l = sp.asarray(sp.shape(self._fibre_image))
        px = sp.zeros(l[0])
        py = sp.zeros(l[1])
        pz = sp.zeros(l[2])

        for x in sp.arange(l[0]):
            px[x] = sp.sum(self._fibre_image[x, :, :])
            px[x] /= sp.size(self._fibre_image[x, :, :])
        for y in sp.arange(l[1]):
            py[y] = sp.sum(self._fibre_image[:, y, :])
            py[y] /= sp.size(self._fibre_image[:, y, :])
        for z in sp.arange(l[2]):
            pz[z] = sp.sum(self._fibre_image[:, :, z])
            pz[z] /= sp.size(self._fibre_image[:, :, z])

        if fig is None:
            fig = plt.figure()
        ax = fig.gca()
        plots = []
        plots.append(plt.plot(sp.arange(l[0])/l[0], px, 'r', label='x'))
        plots.append(plt.plot(sp.arange(l[1])/l[1], py, 'g', label='y'))
        plots.append(plt.plot(sp.arange(l[2])/l[2], pz, 'b', label='z'))
        plt.xlabel('Normalized Distance')
        plt.ylabel('Porosity')
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels, loc=1)
        plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)
        return fig

    def compress_geometry(self, factor=None, preserve_fibres=False):
        r"""
        Adjust the vertices and recalculate geometry. Save fibre voxels before
        and after then put them back into the image to preserve fibre volume.
        Shape will not be conserved. Also make adjustments to the pore and throat
        properties given an approximate volume change from adding the missing fibre
        voxels back in

        Parameters
        ----------
        factor : array_like
        List of 3 values, [x,y,z], 2 must be one and the other must be between
        zero and one representing the fraction of the domain height to compress
        to.

        preserve_fibres : boolean
        If the fibre image has been generated and used to calculate pore volumes
        then preserve fibre volume artificially by adjusting pore and throat sizes
        """
        if factor is None:
            logger.warning('Please supply a compression factor in the form ' +
                           '[1,1,CR], with CR < 1')
            return
        if sp.any(sp.asarray(factor) > 1):
            logger.warning('The supplied compression factor is greater than 1, ' +
                           'the method is not tested for stretching')
            return
        # uncompressed number of fibre voxels in each pore
        fvu = self["pore.fibre_voxels"]
        r1 = self["pore.diameter"]/2
        # Most likely boundary pores - prevents divide by zero (vol change zero)
        r1[r1 == 0] = 1
        vo.scale(network=self._net, scale_factor=factor, preserve_vol=False)
        self.models.regenerate()

        if preserve_fibres and self._voxel_vol:
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
        else:
            logger.warning('Fibre volume is not be conserved under compression')
        # Remove pores with zero throats
        topotools.trim_occluded_throats(self._net)
