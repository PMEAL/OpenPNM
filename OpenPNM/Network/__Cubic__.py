"""
module __Cubic__: Generate lattice-like networks
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry/__init__.py' file.

"""
import OpenPNM

import numpy as np
import scipy as sp
import scipy.sparse as sprs
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
    def __init__(self, shape=None, template=None, spacing=1, **kwargs):
        super(Cubic, self).__init__(**kwargs)
        
        if shape != None:
            arr = np.atleast_3d(np.empty(shape))
        elif template != None:
            arr = sp.array(template,ndmin=3,dtype=bool)
        
        self._shape = sp.shape(arr)  # Store original network shape
        self._spacing = spacing  # Store network spacing instead of calculating it
        
        points = np.array([i for i,v in np.ndenumerate(arr)], dtype=float)
        points += 0.5
        points *= spacing

        I = np.arange(arr.size).reshape(arr.shape)
        tails, heads = [], []
        for T,H in [
            (I[:,:,:-1], I[:,:,1:]),
            (I[:,:-1], I[:,1:]),
            (I[:-1], I[1:]),
            ]:
            tails.extend(T.flat)
            heads.extend(H.flat)
        pairs = np.vstack([tails, heads]).T
        
        self['pore.coords'] = points
        self['throat.conns'] = pairs
        self['pore.all'] = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all'] = np.ones(len(self['throat.conns']), dtype=bool)
        self['pore.index'] = sp.arange(0,len(self['pore.coords']))

        x,y,z = self['pore.coords'].T
        self['pore.internal'] = self['pore.all']
        self['pore.front'] = x <= x.min()
        self['pore.back'] = x >= x.max()
        self['pore.left'] = y <= y.min()
        self['pore.right'] = y >= y.max()
        self['pore.bottom'] = z <= z.min()
        self['pore.top'] = z >= z.max()
        
        #If an image was sent as 'template', then trim network to image shape
        if template != None:
            self.trim(~arr.flatten())

    def add_boundaries(self):
        r'''
        This method uses ``clone`` to clone the surface pores (labeled 'left',
        'right', etc), then shifts them to the periphery of the domain, and
        gives them the label 'right_face', 'left_face', etc.
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
            self.clone(pores=ps,apply_label=[label+'_boundary','boundary'])
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
            self._logger.warning('The supplied pores are not coplanar. Area will be approximate')
        return A
    
    
    
    