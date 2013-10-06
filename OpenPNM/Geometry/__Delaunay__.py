"""
module __Delaunay__: Generate random networks based on Delaunay Tessellations
==========================================================

.. warning:: The classes of this module should be loaded through the 'Geometry.__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.ndimage as spim

from __GenericGeometry__ import GenericGeometry

class Delaunay(GenericGeometry):
    r"""
    Delaunay - Class to create a random network based on the Delaunay tessellation
    
    Parameters
    ----------
    
    domain_size : list with 3 float elements 
        Shape of the domain [Dx, Dy, Dz]
    num_pores :  integer
        Number of randomly located pores
    loglevel : int
        Level of the logger (10=Debug, 20=Info, 30=Warning, 40=Error, 50=Critical)
        
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
        
    """
    
    def __init__(self,**kwargs):
        r"""
        
        """
        super(Delaunay,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        #Instantiate the network
        self._net=OpenPNM.Network.GenericNetwork()

    def _generate_setup(self, **params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")
        self.domain_size = params['domain_size']
        self.num_pores = params['num_pores']
    
    def _generate_pores(self):
        r"""
        Generate the pores with numbering scheme.
        """
        self._logger.info("generate_pores: Place randomly located pores in the domain")
        Lx = self.domain_size[0]
        Ly = self.domain_size[1]
        Lz = self.domain_size[2]
        Np = self.num_pores
        self._net.pore_properties['coords'] = np.hstack((np.random.rand(Np,1)*Lx, np.random.rand(Np,1)*Ly, np.random.rand(Np,1)*Lz))
        self._net.pore_properties['numbering'] = np.arange(0,Np)
        self._net.pore_properties['type']= np.zeros((Np,))
        self._logger.debug("generate_pores: End of method")
        
    def _generate_throats(self):
        r"""
        Generate the throats (connections, numbering and types)
        """
        self._logger.info("generate_throats: Define connections between pores")
        Np = self._net.get_num_pores()
        pts = self._net.pore_properties['coords']
        #Generate 6 dummy domains to pad onto each face of real domain
        #This prevents surface pores from making long range connections to each other
        Lx = self.domain_size[0]
        Ly = self.domain_size[1]
        Lz = self.domain_size[2]
        f = 0.1; #Scale factor for size of dummy domains
        ptsX0 = np.hstack((-np.random.rand(Np*f,1)*Lx*f, np.random.rand(Np*f,1)*Ly, np.random.rand(Np*f,1)*Lz))
        ptsY0 = np.hstack((np.random.rand(Np*f,1)*Lx, -np.random.rand(Np*f,1)*Ly*f, np.random.rand(Np*f,1)*Lz))
        ptsZ0 = np.hstack((np.random.rand(Np*f,1)*Lx, np.random.rand(Np*f,1)*Ly, -np.random.rand(Np*f,1)*Lz*f))
        ptsXX = np.hstack((np.random.rand(Np*f,1)*Lx*f+Lx, np.random.rand(Np*f,1)*Ly, np.random.rand(Np*f,1)*Lz))
        ptsYY = np.hstack((np.random.rand(Np*f,1)*Lx, np.random.rand(Np*f,1)*Ly*f+Ly, np.random.rand(Np*f,1)*Lz))
        ptsZZ = np.hstack((np.random.rand(Np*f,1)*Lx, np.random.rand(Np*f,1)*Ly, np.random.rand(Np*f,1)*Lz*f+Lz))
        #Add dummy domains to real domain
        pts = np.concatenate([pts,ptsX0,ptsXX,ptsY0,ptsYY,ptsZ0,ptsZZ])
        #Perform tessellation
        Tri = sptl.Delaunay(pts)
        adjmat = sprs.lil_matrix((Np,Np),dtype=int)
        for i in np.arange(0,np.shape(Tri.simplices)[0]):
            #Keep only simplices that are fully in real domain
            adjmat[Tri.simplices[i][Tri.simplices[i]<Np],Tri.simplices[i][Tri.simplices[i]<Np]] = 1
        #Remove duplicate (lower triangle) and self connections (diagonal) 
        #and convert to coo
        adjmat = sprs.triu(adjmat,k=1,format="coo")
        self._net.throat_properties['connections'] = np.vstack((adjmat.row, adjmat.col)).T
        self._net.throat_properties['type'] = np.zeros_like(adjmat.row)
        self._net.throat_properties['numbering'] = np.arange(0,np.size(adjmat.row))
        self._logger.debug("generate_throats: End of method")
        
    def _add_boundaries(self):
        r"""
        This is an alternative means of adding boundaries
        """
        self._logger.info("add_boundaries: start of method")
        
        import scipy.spatial as sptl
        import scipy.sparse as sprs
        Lx = self.domain_size[0]
        Ly = self.domain_size[1]
        Lz = self.domain_size[2]
        Np = self._net.get_num_pores()
        btype = self._btype
        boffset = 0.05
        
        #Translate internal pores to each face of domain
        poffset = np.zeros((7,3))
        poffset[[2,5],0] = [-Lx, Lx]
        poffset[[3,4],1] = [-Ly, Ly]
        poffset[[1,6],2] = [-Lz, Lz]
        pcoords = pcoords0 = self._net.pore_properties['coords']
        for i in np.r_[1:7]:
            pcoords = np.concatenate((pcoords,pcoords0 + poffset[i,:]),axis=0)
        
        #Use some twisted logic to get bval list of + for boundary and - for periodic faces
        bval = [0, 1, 2, 3, 4, 5, 6]*(np.array([0, btype[2], btype[0], btype[1], btype[1], btype[0], btype[2]])*-2+1)
        ptype = np.zeros((Np,),dtype=int)
        for i in np.r_[1:7]:
            ptype = np.concatenate((ptype,np.ones((Np,),dtype=int)*bval[i]),axis=0)
        
        #pnum contains the internal ID number of the boundary pores (for connecting periodic points)
        pnum = self._net.pore_properties['numbering'][self._net.pore_properties['type']==0]
        pnum = np.tile(pnum,7)
        
        Tri = sptl.Delaunay(pcoords)
        adjmat = sprs.lil_matrix((np.shape(pcoords)[0],np.shape(pcoords)[0]),dtype=int)
        for i in np.arange(0,np.shape(Tri.simplices)[0]):
            #Keep only simplices that are fully in real domain
            adjmat[Tri.simplices[i],Tri.simplices[i]] = 1
        adjmat = sprs.triu(adjmat,k=1,format="lil")
        for i in np.arange(0,Np):
            #Add periodic throats to the netowrk (if any)
            tpore2 = pnum[adjmat.rows[i]][ptype[adjmat.rows[i]]<0]
            tpore1 = np.ones_like(tpore2,dtype=int)*i
            self._net.throat_properties['connections'] = np.concatenate((self._net.throat_properties['connections'],np.vstack((tpore1,tpore2)).T),axis=0)
            self._net.throat_properties['type'] = np.concatenate((self._net.throat_properties['type'],ptype[adjmat.rows[i]][ptype[adjmat.rows[i]]<0]))
            #Add boundary pores and throats to the network
            newporetyps = np.unique(ptype[adjmat.rows[i]][ptype[adjmat.rows[i]]>0])
            newporenums = np.r_[self._net.get_num_pores():self._net.get_num_pores()+np.size(newporetyps)]
            tpore2 = newporenums
            tpore1 = np.ones_like(tpore2,dtype=int)*i
            self._net.throat_properties['connections'] = np.concatenate((self._net.throat_properties['connections'],np.vstack((tpore1,tpore2)).T),axis=0)
            self._net.throat_properties['type'] = np.concatenate((self._net.throat_properties['type'],newporetyps),axis=0)            
            self._net.pore_properties['type'] = np.concatenate((self._net.pore_properties['type'],newporetyps),axis=0)
            bcoords = np.zeros((7,3),dtype=float)
            bcoords[1,:] = [self._net.pore_properties['coords'][i,0], self._net.pore_properties['coords'][i,1], 0-Lz*boffset]
            bcoords[2,:] = [0-Lx*boffset, self._net.pore_properties['coords'][i,1], self._net.pore_properties['coords'][i,2]]
            bcoords[3,:] = [self._net.pore_properties['coords'][i,0], -Ly*boffset, self._net.pore_properties['coords'][i,2]]
            bcoords[4,:] = [self._net.pore_properties['coords'][i,0], Ly+Ly*boffset, self._net.pore_properties['coords'][i,2]]
            bcoords[5,:] = [Lx+Lx*boffset, self._net.pore_properties['coords'][i,1], self._net.pore_properties['coords'][i,2]]
            bcoords[6,:] = [self._net.pore_properties['coords'][i,0], self._net.pore_properties['coords'][i,1], Lz+Lz*boffset]
            newporecoords = bcoords[newporetyps,:]
            self._net.pore_properties['coords'] = np.concatenate((self._net.pore_properties['coords'],newporecoords),axis=0)
        #Reset number of pores and throats (easier than tracking it)
        self._net.pore_properties['numbering'] = np.r_[0:np.size(self._net.pore_properties['type'])]
        self._net.throat_properties['numbering'] = np.r_[0:np.size(self._net.throat_properties['type'])]
        self._logger.debug("add_boundaries: end of method")
    
    def _add_boundaries_old(self):
        self._logger.info("add_boundaries_old: Start of method")
        
        self.add_opposing_boundaries(btype=[2,5])
        self.add_opposing_boundaries(btype=[3,4])
        self.add_opposing_boundaries(btype=[1,6])    
        
    def _add_opposing_boundaries(self,btype=[1,6]):
        r"""
        btype indicates which two boundaries are being added by type
        """
        self._logger.info("add_opposing_boundaries: start of method") 
        
        if btype==[2,5]:
            D=0
            W=1
            H=2
        elif btype==[3,4]:
            D=1
            W=0
            H=2
        elif btype==[1,6]:
            D=2
            W=1
            H=0

        Lx = self.domain_size[D]
        Ly = self.domain_size[W]
        Lz = self.domain_size[H]
        #Rotate pore coordinates (use only internal pores)
        pnum = self._net.pore_properties['numbering'][self._net.pore_properties['type']==0]
        pcoords = np.zeros_like(self._net.pore_properties['coords'][pnum,:])
        pcoords[:,0] = self._net.pore_properties['coords'][pnum,D]
        pcoords[:,1] = self._net.pore_properties['coords'][pnum,W]
        pcoords[:,2] = self._net.pore_properties['coords'][pnum,H]
        
        #Determine dimensions of image from dimensions of domain
        f = 100 #minimum image dimension
        im_dim = [0,0,0]
        im_dim[0] = np.floor(f*Lx/np.min([Lx,Ly,Lz]))
        im_dim[1] = np.floor(f*Ly/np.min([Lx,Ly,Lz]))
        im_dim[2] = np.floor(f*Lz/np.min([Lx,Ly,Lz]))
        im_dim = np.array(im_dim,dtype=int)
        
        #Convert pore coordinates into image subscripts
        im_subs = np.zeros_like(pcoords,dtype=int)
        im_subs[:,0] = pcoords[:,0]*im_dim[0]/Lx
        im_subs[:,1] = pcoords[:,1]*im_dim[1]/Ly
        im_subs[:,2] = pcoords[:,2]*im_dim[2]/Lz
        #Find linear indices of each pore in the new image
        im_inds = np.ravel_multi_index((im_subs[:,0], im_subs[:,1], im_subs[:,2]), dims=(im_dim), order='F')
        
        #Generate 3D image of points (place pore numbers at each site for use later)
        img = np.zeros(im_dim,dtype=int)
        img[im_subs[:,0],im_subs[:,1],im_subs[:,2]] = pnum

        #Perform distance transform on points and also get 'indicies' of each point
        img_dt, ind_dt = spim.distance_transform_edt(img==0,return_indices=True)
        
        #Project all* internal points to x face
        #*Note that it's possible/likely that mutliple internal points map to the same boundary point
        img_bd0 = np.zeros([im_dim[1],im_dim[2]],dtype=int)
        img_bd1 = np.zeros([im_dim[1],im_dim[2]],dtype=int)
        img_bd0[im_subs[:,1],im_subs[:,2]] = im_inds
        img_bd1[im_subs[:,1],im_subs[:,2]] = im_inds
        
        #Create 2D array of distance transform indices for 0 and end faces
        dt_D0 = ind_dt[0,0,:,:]*(img_bd0>0)  #0 face
        dt_D1 = ind_dt[0,-1,:,:]*(img_bd1>0) #end face
        
        #Create a 2D mask containing x coordinates of internal points and -1 elsewhere
        img_D0 = -np.ones([im_dim[1],im_dim[2]],dtype=int)
        img_D1 = -np.ones([im_dim[1],im_dim[2]],dtype=int)
        img_D0[im_subs[:,1],im_subs[:,2]] = im_subs[:,0]
        img_D1[im_subs[:,1],im_subs[:,2]] = im_subs[:,0]
        
        #Find where x value of internal points corresponds to x value of distance transform indices
        img_bd0 = (img_D0 == dt_D0)*img_bd0
        img_bd1 = (img_D1 == dt_D1)*img_bd1
        
        #Convert boundary sites to linear indices
        inds_bd0 = img_bd0[np.nonzero(img_bd0)]
        inds_bd1 = img_bd1[np.nonzero(img_bd1)]
        
        #Use linear indices to find pore ID nums
        nums_bd0 = img[np.unravel_index(inds_bd0, dims=(im_dim), order='F')]
        nums_bd1 = img[np.unravel_index(inds_bd1, dims=(im_dim), order='F')]
        nums_bd = np.append(nums_bd0,nums_bd1)
        types_bd = np.append(np.zeros_like(nums_bd0),np.ones_like(nums_bd1))
        
        #Add new boundary pores and throats to the network
        Np = self._net.get_num_pores() #Get all pores including previously added boundaries
        bp_numbering = np.r_[Np:Np+np.size(nums_bd)]
        bp_type = (types_bd==0)*btype[0] + (types_bd==1)*btype[1]
        bp_coords = np.zeros([np.size(nums_bd),3])
        bp_coords[types_bd==0,D] = np.zeros_like(nums_bd0)-.0001
        bp_coords[types_bd==0,W] = pcoords[nums_bd0,1]
        bp_coords[types_bd==0,H] = pcoords[nums_bd0,2]
        bp_coords[types_bd==1,D] = np.ones_like(nums_bd1)*Lx+0.0001
        bp_coords[types_bd==1,W] = pcoords[nums_bd1,1]
        bp_coords[types_bd==1,H] = pcoords[nums_bd1,2]
        self._net.pore_properties['numbering'] = np.append(self._net.pore_properties['numbering'],bp_numbering)
        self._net.pore_properties['type'] = np.append(self._net.pore_properties['type'],bp_type)
        self._net.pore_properties['coords'] = np.concatenate((self._net.pore_properties['coords'],bp_coords))
        Nt = self._net.get_num_throats()
        bt_numbering = np.r_[Nt:Nt+np.size(nums_bd)]
        bt_type = np.ones(np.size(nums_bd),dtype=int)*2
        bt_connections = np.zeros([np.size(nums_bd),2],dtype=int)
        bt_connections[:,0] = nums_bd
        bt_connections[:,1] = bp_numbering
        self._net.throat_properties['numbering'] = np.append(self._net.throat_properties['numbering'],bt_numbering)
        self._net.throat_properties['type'] = np.append(self._net.throat_properties['type'],bt_type)
        self._net.throat_properties['connections'] =  np.concatenate((self._net.throat_properties['connections'],bt_connections))
        
if __name__ == '__main__':
    test=Delaunay(loggername='TestDelaunay')
#    test.generate(num_pores=250, domain_size=[100,100,100])