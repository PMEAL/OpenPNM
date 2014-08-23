"""
module __Cubic__: Generate lattice-like networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry/__init__.py' file.

"""
import OpenPNM

import numpy as np
import scipy as sp
import OpenPNM.Utilities.misc as misc
from OpenPNM.Network.__GenericNetwork__ import GenericNetwork

class Cubic(GenericNetwork):
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
    This class is a work in progress, examples forthcoming.
    """
    @classmethod
    def from_image(cls, image, *args, **kwargs):
        network = cls(image.shape, *args, **kwargs)
        network['pore.values'] = image.ravel()
        return network
    
    def __init__(self, shape, spacing=None, bbox=None, **kwargs):
        super(Cubic, self).__init__(**kwargs)
        arr = np.atleast_3d(np.empty(shape))

        points = np.array([i for i,v in np.ndenumerate(arr)], dtype=float)
        points += 0.5
        if spacing is not None:
            points *= spacing
        elif bbox is not None:
            points *= bbox / self.points.max(axis=0)

        I = np.arange(arr.size).reshape(arr.shape)
        tails, heads = [], []
        for T,H in [
            (I[:,:,:-1], I[:,:,1:]),
            (I[:,:-1], I[:,1:]),
            (I[:-1], I[1:]),
            ]:
            tails.extend(T.flat)
            tails.extend(H.flat)
            heads.extend(H.flat)
            heads.extend(T.flat)
        pairs = np.vstack([tails, heads]).T

        self['pore.coords'] = points
        self['throat.conns'] = pairs

        self['pore.all'] = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all'] = np.ones(len(self['throat.conns']), dtype=bool)

        x,y,z = self['pore.coords'].T
        self['pore.internal'] = self['pore.all']
        self['pore.front'] = x <= x.min()
        self['pore.back'] = x >= x.max()
        self['pore.left'] = y <= y.min()
        self['pore.right'] = y >= y.max()
        self['pore.bottom'] = z <= z.min()
        self['pore.top'] = z >= z.max()

    def add_boundaries(self):
        r'''
        This method uses ``clone`` to clone the surface pores (labeled 'left',
        'right', etc), then shifts them to the periphery of the domain, and
        gives them the label 'right_face', 'left_face', etc.
        '''
        x,y,z = self['pore.coords'].T
        
        Lc = sp.amax(sp.diff(x))
        
        offset = {}
        offset['front'] = offset['left'] = offset['bottom'] = [0,0,0]
        offset['back']  = [x.max()+Lc/2,0,0]
        offset['right'] = [0,y.max()+Lc/2,0]
        offset['top']   = [0,0,z.max()+Lc/2]
        
        scale = {}
        scale['front']  = scale['back']  = [0,1,1]
        scale['left']   = scale['right'] = [1,0,1]
        scale['bottom'] = scale['top']   = [1,1,0]
        
        for label in ['front','back','left','right','bottom','top']:
            ps = self.pores(label)
            self.clone(pores=ps,apply_label=[label+'_boundary','boundary'])
            #Translate cloned pores
            ind = self.pores(label+'_boundary')
            coords = self['pore.coords'][ind]
            coords = coords*scale[label] + offset[label]
            self['pore.coords'][ind] = coords

    def asarray(self, values):
        points = self['pore.coords']
        spacing = map(np.diff, map(np.unique, points.T))
        min_spacing = [min(a) if len(a) else 1.0 for a in spacing]
        points = (points / min_spacing).astype(int)
        bbox = points.max(axis=0) - points.min(axis=0)
        bbox = (bbox / min_spacing + 1).astype(int)
        actual_indexes = np.ravel_multi_index(points.T, bbox)
        print(bbox)
        array = np.zeros(bbox)
        array.flat[actual_indexes] = values.ravel()
        return array.squeeze()
        
    def domain_length(self,face_1,face_2):
        r'''
        Calculate the distance between two faces of the network
        
        Parameters
        ----------
        face_1 and face_2 : array_like
            Lists of pores belonging to opposite faces of the network
            
        Returns
        -------
        The length of the domain in the specified direction
        
        Notes
        -----
        - Does not yet check if input faces are perpendicular to each other
        '''
        #Ensure given points are coplanar before proceeding
        if misc.iscoplanar(self['pore.coords'][face_1]) and misc.iscoplanar(self['pore.coords'][face_2]):
            #Find distance between given faces
            x = self['pore.coords'][face_1]
            y = self['pore.coords'][face_2]
            Ds = misc.dist(x,y)
            L = sp.median(sp.amin(Ds,axis=0))
        else:
            self._logger.warning('The supplied pores are not coplanar. Length will be approximate.')
            f1 = self['pore.coords'][face_1]
            f2 = self['pore.coords'][face_2]
            distavg = [0,0,0]
            distavg[0] = sp.absolute(sp.average(f1[:,0]) - sp.average(f2[:,0]))
            distavg[1] = sp.absolute(sp.average(f1[:,1]) - sp.average(f2[:,1]))
            distavg[2] = sp.absolute(sp.average(f1[:,2]) - sp.average(f2[:,2]))
            L = max(distavg)
        return L

        
    def domain_area(self,face):
        r'''
        Calculate the area of a given network face
        
        Parameters
        ----------
        face : array_like
            List of pores of pore defining the face of interest
            
        Returns
        -------
        The area of the specified face
        '''
        coords = self['pore.coords'][face]
        dx = max(coords[:,0]) - min(coords[:,0])
        dy = max(coords[:,1]) - min(coords[:,1])
        dz = max(coords[:,2]) - min(coords[:,2])
        xy = dx*dy
        xz = dx*dz
        yz = dy*dz
        A = max([xy,xz,yz])
        
        if not misc.iscoplanar(self['pore.coords'][face]):
            self._logger.warning('The supplied pores are not coplanar. Area will be approximate')
        return A
    
    
    
    