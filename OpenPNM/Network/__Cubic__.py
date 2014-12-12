# -*- coding: utf-8 -*-
"""
===============================================================================
Cubic: Generate lattice-like networks
===============================================================================

"""
import numpy as np
import scipy as sp
import OpenPNM.Utilities.misc as misc
from OpenPNM.Network import GenericNetwork
from OpenPNM.Base import logging
logger = logging.getLogger()

class Cubic(GenericNetwork):
    r"""
    This class generates a cubic network of the specified size and shape.
    Alternatively, an arbitrary domain shape defined by a supplied template.

    Parameters
    ----------
    name : string
        A unique name for the network

    shape : tuple of ints
        The (i,j,k) size and shape of the network.

    template : array of booleans
        An (i,j,k) array with True where the Network should be defiend and
        False elsewhere.  This approach is useful for creating networks of non-
        cuboid shape like spheres or cylinders, but still with a cubic lattice
        topology.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.Cubic(shape=[3,4,5])
    >>> pn.Np
    60
    >>> img = sp.ones([3,4,5])
    >>> pn = OpenPNM.Network.Cubic(template=img)
    >>> pn.Np
    60
    """
    def __init__(self, shape=None, template=None, spacing=1, connectivity=6, **kwargs):
        super(Cubic, self).__init__(**kwargs)

        if shape is not None:
            arr = np.atleast_3d(np.empty(shape))
        elif template is not None:
            arr = sp.array(template,ndmin=3,dtype=bool)
        else:
            arr = np.atleast_3d(np.empty([1,1,1]))

        self._shape = sp.shape(arr)  # Store original network shape
        self._spacing = spacing  # Store network spacing instead of calculating it

        points = np.array([i for i,v in np.ndenumerate(arr)], dtype=float)
        points += 0.5
        points *= spacing

        I = np.arange(arr.size).reshape(arr.shape)

        face_joints = [
            (I[:,:,:-1], I[:,:,1:]),
            (I[:,:-1], I[:,1:]),
            (I[:-1], I[1:]),
            ]

        edge_joints = [
            (I[:-1,:-1,:-1], I[1:,1:,1:]),
            (I[:-1,:-1,1:], I[1:,1:,:-1]),
            (I[:-1,1:,:-1], I[1:,:-1,1:]),
            (I[1:,:-1,:-1], I[:-1,1:,1:]),
            ]

        corner_joints = [
            (I[:,:-1,:-1], I[:,1:,1:]),
            (I[:,:-1,1:], I[:,1:,:-1]),
            (I[:-1,:,:-1], I[1:,:,1:]),
            (I[1:,:,:-1], I[:-1,:,1:]),
            (I[1:,1:,:], I[:-1,:-1,:]),
            (I[1:,:-1,:], I[:-1,1:,:]),
            ]
        
        if connectivity == 6:
            joints = face_joints
        elif connectivity == 8:
            joints = edge_joints
        elif connectivity == 18:
            joints = face_joints + edge_joints
        elif connectivity == 26:
            joints = face_joints + corner_joints + edge_joints
        else:
            raise Exception('Invalid connectivity receieved. Must be 6, 8, 18 or 26')

        I = np.arange(arr.size).reshape(arr.shape)
        tails, heads = [], []
        for T,H in joints:
            tails.extend(T.flat)
            heads.extend(H.flat)
        pairs = np.vstack([tails, heads]).T

        self['pore.coords']  = points
        self['throat.conns'] = pairs
        self['pore.all']     = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all']   = np.ones(len(self['throat.conns']), dtype=bool)
        self['pore.index']   = sp.arange(0,len(self['pore.coords']))

        x,y,z = self['pore.coords'].T
        self['pore.internal'] = self['pore.all']
        self['pore.front']    = x <= x.min()
        self['pore.back']     = x >= x.max()
        self['pore.left']     = y <= y.min()
        self['pore.right']    = y >= y.max()
        self['pore.bottom']   = z <= z.min()
        self['pore.top']      = z >= z.max()

        #Add some topology models to the Network
#        mod = OpenPNM.Network.models.pore_topology.get_subscripts
#        self.add_model(propname='pore.subscript',model=mod,shape=self._shape)

        #If an image was sent as 'template', then trim network to image shape
        if template is not None:
            self.trim(~arr.flatten())

    def add_boundaries(self):
        r'''
        This method uses ``clone_pores`` to clone the surface pores (labeled
        'left','right', etc), then shifts them to the periphery of the domain,
        and gives them the label 'right_face', 'left_face', etc.
        '''
        x,y,z = self['pore.coords'].T

        Lc = sp.amax(sp.diff(x)) #this currently works but is very fragile

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
            self.clone_pores(pores=ps,apply_label=[label+'_boundary','boundary'])
            #Translate cloned pores
            ind = self.pores(label+'_boundary')
            coords = self['pore.coords'][ind]
            coords = coords*scale[label] + offset[label]
            self['pore.coords'][ind] = coords

    def asarray(self,values):
        r'''
        Retreive values as a rectangular array, rather than the OpenPNM list format

        Parameters
        ----------
        values : array_like
            The values from the network (in a list) to insert into the array

        Notes
        -----
        This method can break on networks that have had boundaries added.  It
        will usually work IF the list of values came only from 'internal' pores.
        '''
        if sp.shape(values)[0] > self.num_pores('internal'):
            raise Exception('The received values are bigger than the original network')
        Ps = sp.array(self['pore.index'][self.pores('internal')],dtype=int)
        arr = sp.ones(self._shape)*sp.nan
        ind = sp.unravel_index(Ps,self._shape)
        arr[ind[0],ind[1],ind[2]] = values
        return arr

    def fromarray(self,array,propname):
        r'''
        Apply data to the network based on a rectangular array filled with
        values.  Each array location corresponds to a pore in the network.

        Parameters
        ----------
        array : array_like
            The rectangular array containing the values to be added to the
            network. This array must be the same shape as the original network.

        propname : string
            The name of the pore property being added.
        '''
        array = sp.atleast_3d(array)
        if sp.shape(array) != self._shape:
            raise Exception('The received array does not match the original network')
        temp = array.flatten()
        Ps = sp.array(self['pore.index'][self.pores('internal')],dtype=int)
        propname = 'pore.' + propname.split('.')[-1]
        self[propname] = sp.nan
        self[propname][self.pores('internal')] = temp[Ps]

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
            logger.warning('The supplied pores are not coplanar. Length will be approximate.')
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
        rads = self['pore.diameter'][face]/2.
        # calculate the area of the 3 principle faces of the bounding cuboid
        dx = max(coords[:,0]+rads) - min(coords[:,0]-rads)
        dy = max(coords[:,1]+rads) - min(coords[:,1]-rads)
        dz = max(coords[:,2]+rads) - min(coords[:,2]-rads)
        yz = dy*dz # x normal
        xz = dx*dz # y normal
        xy = dx*dy # z normal
        # find the directions parallel to the plane
        directions = sp.where([yz,xz,xy]!=max([yz,xz,xy]))[0]
        try:
            # now, use the whole network to do the area calculation
            coords = self['pore.coords']
            rads = self['pore.diameter']/2.
            d0 = (max(coords[:,directions[0]]+rads) - min(coords[:,directions[0]]-rads))
            d1 = (max(coords[:,directions[1]]+rads) - min(coords[:,directions[1]]-rads))
            A = d0*d1
        except:
            # if that fails, use the max face area of the bounding cuboid
            A = max([yz,xz,xy])
        if not misc.iscoplanar(self['pore.coords'][face]):
            logger.warning('The supplied pores are not coplanar. Area will be approximate')
            pass
        return A


if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)

