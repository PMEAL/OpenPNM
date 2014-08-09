"""
module __Delaunay__: Generate random networks based on Delaunay Tessellations
==========================================================

.. warning:: The classes of this module should be loaded through the 'Topology.__init__.py' file.

"""

import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM

import sys
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.ndimage as spim
from OpenPNM.Network.__GenericNetwork__ import GenericNetwork
from scipy.spatial import Voronoi

class Delaunay(GenericNetwork):
    r"""
    This class contains the methods for creating a *Delaunay* network topology
    based connecting pores with a Delaunay tessellation.  
    
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
    >>> pn = OpenPNM.Network.Delaunay()
    >>> pn.generate(num_pores=100,domain_size=[100,100,100],add_boundaries=0)
    >>> pn.num_pores()
    100

    """
    add_boundaries = False

    def __init__(self,**kwargs):
        '''
        Create Delauny network object
        '''
        super(Delaunay,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
        
    def generate(self,**params):
        r'''
        Method to trigger the generation of the network
        
        Parameters
        ----------
        domain_size : list of floats, [Lx,Ly,Lz]
            Bounding cube for internal pore positions
        num_pores : int
            Number of pores to place randomly within domain

        '''
        self._logger.info(sys._getframe().f_code.co_name+": Start of network topology generation")
        self._generate_setup(**params)
        if params['add_boundaries']:
            self.add_boundaries = True
        self._generate_pores()
        self._generate_throats()
        #self._add_labels()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")

    def _generate_setup(self, **params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")
        if params['domain_size'] and params['num_pores']:
            self._Lx = params['domain_size'][0]
            self._Ly = params['domain_size'][1]
            self._Lz = params['domain_size'][2]
            self._Np = params['num_pores']
            r'''
            TODO: Fix this, btype should be received as an argument
            '''
            self._btype = [0,0,0]
        else:
            self._logger.error("domain_size and num_pores must be specified")
            raise Exception('domain_size and num_pores must be specified')

    def _generate_pores(self):
        r"""
        Generate the pores with numbering scheme.
        """
        self._logger.info(sys._getframe().f_code.co_name+": Place randomly located pores in the domain")
        coords = sp.rand(self._Np,3)*[self._Lx,self._Ly,self._Lz]
        self.set_pore_data(prop='coords',data=coords)
        self.set_pore_info(label='all',locations=np.ones_like(coords[:,0]))
        self._logger.debug(sys._getframe().f_code.co_name+": End of method")    
    
    def _generate_throats(self):
        r"""
        Generate the throats connections
        """
        self._logger.info(sys._getframe().f_code.co_name+": Define connections between pores")
        Np = self.num_pores()
        pts = self.get_pore_data(prop='coords')
        #Generate 6 dummy domains to pad onto each face of real domain
        #This prevents surface pores from making long range connections to each other
        Lx = self._Lx
        Ly = self._Ly
        Lz = self._Lz
        #f = 0.1; #Scale factor to reduce size of dummy domains
        #Np_f = sp.array(Np*f,dtype=int)
        #ptsX0 = sp.rand(Np_f,3)*sp.array([-Lx*f,Ly*f,Lz*f])
        #ptsY0 = sp.rand(Np_f,3)*[Lx*f,-Ly*f,Lz*f]
        #ptsZ0 = sp.rand(Np_f,3)*[Lx*f,Ly*f,-Lz*f]
        #ptsXX = sp.rand(Np_f,3)*[Lx*f,Ly*f,Lz*f]+[Lx,0,0]
        #ptsYY = sp.rand(Np_f,3)*[Lx*f,Ly*f,Lz*f]+[0,Ly,0]
        #ptsZZ = sp.rand(Np_f,3)*[Lx*f,Ly*f,Lz*f]+[0,0,Lz]
        
        " Reflect in X = Lx and 0 "
        Pxp = pts.copy()
        Pxp[:,0]=(2*Lx-Pxp[:,0])
        Pxm= pts.copy()
        Pxm[:,0] = Pxm[:,0]*(-1)
        " Reflect in Y = Ly and 0 "
        Pyp = pts.copy()
        Pyp[:,1]=(2*Ly-Pxp[:,1])
        Pym = pts.copy()
        Pym[:,1] = Pxm[:,1]*(-1)
        " Reflect in Z = Lz and 0 "
        Pzp = pts.copy()
        Pzp[:,2]=(2*Lz-Pxp[:,2])
        Pzm = pts.copy()
        Pzm[:,2] = Pxm[:,2]*(-1)
        pts = np.vstack((pts,Pxp,Pxm,Pyp,Pym,Pzp,Pzm)) #Order important for boundary logic
        #Add dummy domains to real domain
        #pts = sp.concatenate([pts,ptsX0,ptsXX,ptsY0,ptsYY,ptsZ0,ptsZZ])
        #Perform tessellation
        self._logger.debug(sys._getframe().f_code.co_name+": Beginning tessellation")
        Tri = sptl.Delaunay(pts)
        self._logger.debug(sys._getframe().f_code.co_name+": Converting tessellation to adjacency matrix")
        adjmat = sprs.lil_matrix((Np,Np),dtype=int)
        for i in sp.arange(0,sp.shape(Tri.simplices)[0]):
            #Keep only simplices that are fully in real domain
            #this used to be vectorize, but it stopped working...change in scipy?
            for j in Tri.simplices[i]: 
                if j < Np:
                    adjmat[j,Tri.simplices[i][Tri.simplices[i]<Np]] = 1
        #Remove duplicate (lower triangle) and self connections (diagonal)
        #and convert to coo
        adjmat = sprs.triu(adjmat,k=1,format="coo")
        self._logger.debug(sys._getframe().f_code.co_name+": Conversion to adjacency matrix complete")
        self.set_throat_data(prop='conns',data=sp.vstack((adjmat.row, adjmat.col)).T)
        tpore1 = self.get_throat_data(prop='conns')[:,0]
        self.set_throat_info(label='all',locations=np.ones_like(tpore1))
        
        if self.add_boundaries:
            " New code to identify boundary pores - those that connect to pores inside and outside original set of pores "
            boundary_pores=[False]*Np
            boundary_pore_list = []
            xp = []
            xm = []
            yp = []
            ym = []
            zp = []
            zm = []
            for i in sp.arange(0,sp.shape(Tri.simplices)[0]):
                pores_in = Tri.simplices[i] < Np # Pores in the original domain
                if (sum(pores_in) >= 1) and (sum(pores_in) < len(pores_in)):
                    # We have connections between at least one pore in and out of the domain
                    # Identify which boundary we are connected to. At corners could be more than one
                    pores_xp = (Tri.simplices[i] >= Np) & (Tri.simplices[i] < 2*Np) # Pores in positive x dummy reflection
                    pores_xm = (Tri.simplices[i] >= 2*Np) & (Tri.simplices[i] < 3*Np) # Pores in negative x dummy reflection
                    pores_yp = (Tri.simplices[i] >= 3*Np) & (Tri.simplices[i] < 4*Np) # Pores in positive y dummy reflection
                    pores_ym = (Tri.simplices[i] >= 4*Np) & (Tri.simplices[i] < 5*Np) # Pores in negative y dummy reflection            
                    pores_zp = (Tri.simplices[i] >= 5*Np) & (Tri.simplices[i] < 6*Np) # Pores in positive z dummy reflection
                    pores_zm = (Tri.simplices[i] >= 6*Np) & (Tri.simplices[i] < 7*Np) # Pores in negative z dummy reflection            
                    for j in range(len(Tri.simplices[i])):
                        if pores_in[j] == True:
                            pore_id = Tri.simplices[i][j]
                            boundary_pores[pore_id]=True
                            if pore_id not in boundary_pore_list:
                                boundary_pore_list.append(pore_id)
                            if sum(pores_xp) >= 1 and pore_id not in xp:
                                xp.append(pore_id)
                            if sum(pores_xm) >= 1 and pore_id not in xm:
                                xm.append(pore_id)
                            if sum(pores_yp) >= 1 and pore_id not in yp:
                                yp.append(pore_id)
                            if sum(pores_ym) >= 1 and pore_id not in ym:
                                ym.append(pore_id)
                            if sum(pores_zp) >= 1 and pore_id not in zp:
                                zp.append(pore_id)
                            if sum(pores_zm) >= 1 and pore_id not in zm:
                                zm.append(pore_id)
                            
            vor_bounds=sp.asarray(boundary_pore_list)
            vor_bound_back = sp.asarray(xp)
            vor_bound_front = sp.asarray(xm)
            vor_bound_right = sp.asarray(yp)
            vor_bound_left = sp.asarray(ym)
            vor_bound_top = sp.asarray(zp)
            vor_bound_bottom = sp.asarray(zm)
            self.set_pore_info(label='boundary',locations=vor_bounds)
            self.set_pore_info(label='right',locations=vor_bound_right) 
            self.set_pore_info(label='left',locations=vor_bound_left) 
            self.set_pore_info(label='front',locations=vor_bound_front) 
            self.set_pore_info(label='back',locations=vor_bound_back) 
            self.set_pore_info(label='top',locations=vor_bound_top) 
            self.set_pore_info(label='bottom',locations=vor_bound_bottom)
            self.set_pore_info(label='internal',locations='all')
        
        # Do Voronoi diagram - creating voronoi polyhedra around each pore and save vertex information
        vor = Voronoi(pts)
        all_verts = sp.ndarray(Np,dtype=object)
        for i,polygon in enumerate(vor.point_region[0:Np]):
            if -1 not in vor.regions[polygon]:
                all_verts[i]=vor.vertices[vor.regions[polygon]]
            else:
                all_verts[i]="unbounded"
        self.set_pore_data(prop='vertices',data=all_verts)
        self._logger.debug(sys._getframe().f_code.co_name+": End of method")
        
    def _add_labels(self):
        r'''
        This finds surface pors simply by proximity to the domain boundaries.
        A better approach is necessary 
        '''
        coords = self.get_pore_data(prop='coords')
        self.set_pore_info(label='front',locations=(coords[:,0]<(0.1*self._Lx)))
        self.set_pore_info(label='back',locations=(coords[:,0]>(0.9*self._Lx)))
        self.set_pore_info(label='left',locations=(coords[:,1]<(0.1*self._Ly)))
        self.set_pore_info(label='right',locations=(coords[:,1]>(0.9*self._Ly)))
        self.set_pore_info(label='bottom',locations=(coords[:,2]<(0.1*self._Lz)))
        self.set_pore_info(label='top',locations=(coords[:,2]>(0.9*self._Lz)))
        bnds = self.get_pore_indices(labels=['front','back','left','right','bottom','top'])
        self.set_pore_info(label='boundary',locations=bnds)
        self.set_pore_info(label='internal',locations='all')
        
    def _add_boundaries(self):
        r"""
        This is an alternative means of adding boundaries
        """
        self._logger.info("add_boundaries: start of method")

        import scipy.spatial as sptl
        import scipy.sparse as sprs
        Lx = self._Lx
        Ly = self._Ly
        Lz = self._Lz
        Np = self.num_pores()
        btype = self._btype
        boffset = 0.05

        #Translate internal pores to each face of domain
        poffset = np.zeros((7,3))
        poffset[[2,5],0] = [-Lx, Lx]
        poffset[[3,4],1] = [-Ly, Ly]
        poffset[[1,6],2] = [-Lz, Lz]
        pcoords = pcoords0 = self.get_pore_data(prop='coords')
        for i in np.r_[1:7]:
            pcoords = np.concatenate((pcoords,pcoords0 + poffset[i,:]),axis=0)

        #Use some twisted logic to get bval list of + for boundary and - for periodic faces
        bval = [0, 1, 2, 3, 4, 5, 6]*(np.array([0, btype[2], btype[0], btype[1], btype[1], btype[0], btype[2]])*-2+1)
        ptype = np.zeros((Np,),dtype=int)
        for i in np.r_[1:7]:
            ptype = np.concatenate((ptype,np.ones((Np,),dtype=int)*bval[i]),axis=0)

        #pnum contains the internal ID number of the boundary pores (for connecting periodic points)
        pnum = self.get_pore_indices()
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
            conns = self.get_throat_data(prop='conns')
            conns = np.concatenate((conns,np.vstack((tpore1,tpore2)).T),axis=0)
            #Add boundary pores and throats to the network
            newporetyps = np.unique(ptype[adjmat.rows[i]][ptype[adjmat.rows[i]]>0])
            newporenums = np.r_[self.num_pores():self.num_pores()+np.size(newporetyps)]
            tpore2 = newporenums
            tpore1 = np.ones_like(tpore2,dtype=int)*i
            conns = np.concatenate((conns,np.vstack((tpore1,tpore2)).T),axis=0)
            self.set_throat_data(prop='conns',data=conns)
            bcoords = np.zeros((7,3),dtype=float)
            coords = self.get_pore_data(prop='coords')
            bcoords[1,:] = [coords[i,0], coords[i,1], 0-Lz*boffset]
            bcoords[2,:] = [0-Lx*boffset, coords[i,1], coords[i,2]]
            bcoords[3,:] = [coords[i,0], -Ly*boffset, coords[i,2]]
            bcoords[4,:] = [coords[i,0], Ly+Ly*boffset, coords[i,2]]
            bcoords[5,:] = [Lx+Lx*boffset, coords[i,1], coords[i,2]]
            bcoords[6,:] = [coords[i,0], coords[i,1], Lz+Lz*boffset]
            newporecoords = bcoords[newporetyps,:]
            coords = np.concatenate((coords,newporecoords),axis=0)
            self.set_pore_data(prop='coords',data=coords)
        #Reset number of pores and throats (easier than tracking it)
        nums = np.r_[0:np.shape(coords)[0]]
        self.set_pore_data(prop='numbering',data=nums)
        self.set_pore_info(label='numbering',locations=np.ones((nums[-1]+1,),dtype=bool))
        nums = np.r_[0:np.shape(conns)[0]]
        self.set_throat_data(prop='numbering',data=nums)
        self.set_throat_info(label='numbering',locations=np.ones((nums[-1]+1,),dtype=bool))
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
        pnum = self._net.pore_data['numbering'][self._net.pore_data['type']==0]
        pcoords = np.zeros_like(self._net.pore_data['coords'][pnum,:])
        pcoords[:,0] = self._net.pore_data['coords'][pnum,D]
        pcoords[:,1] = self._net.pore_data['coords'][pnum,W]
        pcoords[:,2] = self._net.pore_data['coords'][pnum,H]

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
        img_dt, ind_dt = spim.distance_transform_edt(img==0)

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
        Np = self._net.num_pores() #Get all pores including previously added boundaries
        bp_numbering = np.r_[Np:Np+np.size(nums_bd)]
        bp_type = (types_bd==0)*btype[0] + (types_bd==1)*btype[1]
        bp_coords = np.zeros([np.size(nums_bd),3])
        bp_coords[types_bd==0,D] = np.zeros_like(nums_bd0)-.0001
        bp_coords[types_bd==0,W] = pcoords[nums_bd0,1]
        bp_coords[types_bd==0,H] = pcoords[nums_bd0,2]
        bp_coords[types_bd==1,D] = np.ones_like(nums_bd1)*Lx+0.0001
        bp_coords[types_bd==1,W] = pcoords[nums_bd1,1]
        bp_coords[types_bd==1,H] = pcoords[nums_bd1,2]
        self._net.pore_data['numbering'] = np.append(self._net.pore_data['numbering'],bp_numbering)
        self._net.pore_data['type'] = np.append(self._net.pore_data['type'],bp_type)
        self._net.pore_data['coords'] = np.concatenate((self._net.pore_data['coords'],bp_coords))
        Nt = self._net.num_throats()
        bt_numbering = np.r_[Nt:Nt+np.size(nums_bd)]
        bt_type = np.ones(np.size(nums_bd),dtype=int)*2
        bt_connections = np.zeros([np.size(nums_bd),2],dtype=int)
        bt_connections[:,0] = nums_bd
        bt_connections[:,1] = bp_numbering
        self._net.throat_data['numbering'] = np.append(self._net.throat_data['numbering'],bt_numbering)
        self._net.throat_data['type'] = np.append(self._net.throat_data['type'],bt_type)
        self._net.throat_data['conns'] =  np.concatenate((self._net.throat_data['conns'],bt_connections))
    
    def domain_size(self,dimension=''):
        if dimension == 'front' or dimension == 'back':
            return self._Ly*self._Lz
        if dimension == 'left' or dimension == 'right':
            return self._Lx*self._Lz
        if dimension == 'top' or dimension == 'bottom':
            return self._Lx*self._Ly
        if dimension == 'volume':
            return self._Lx*self._Ly*self._Lz
        if dimension == 'height':
            return self._Lz
        if dimension == 'width':
            return self._Lx
        if dimension == 'depth':
            return self._Ly
            
if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
