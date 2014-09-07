"""
module __Delaunay__: Generate random networks based on Delaunay Tessellations
==========================================================

.. warning:: The classes of this module should be loaded through the 'Topology.__init__.py' file.

"""
import OpenPNM

import sys
import scipy as sp
import numpy as np
import scipy.sparse as sprs
import scipy.spatial as sptl
import scipy.ndimage as spim
from OpenPNM.Network.__GenericNetwork__ import GenericNetwork
from scipy.spatial import Voronoi
from scipy import stats as st
from scipy.special import cbrt

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
        #self._generate_cubic_pores() - Impose a quasi-cubic pore distribution - needs work
        self._generate_throats()
        #self._add_labels()
        self._logger.debug(sys._getframe().f_code.co_name+": Network generation complete")

    def _generate_setup(self, **params):
        r"""
        Perform applicable preliminary checks and calculations required for generation
        """
        self._logger.debug("generate_setup: Perform preliminary calculations")
        if 'aniso' in params:
            self._aniso = sp.asarray(params['aniso'])
        else:
            self._aniso = sp.array([1,1,1])
        #self._rescale_factor = 1/cbrt(np.prod(self._aniso))
        self._rescale_factor = 1/(self._aniso)
        self._scale_factor = self._aniso*self._rescale_factor
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
        #Original Random Point Generator
        #coords = sp.rand(self._Np,3)*[self._Lx,self._Ly,self._Lz]
        #Introduce Anisotropy
        [self._Lx,self._Ly,self._Lz]=np.around([self._Lx,self._Ly,self._Lz]*self._aniso,10)
        #Seeding Code
        coords = np.zeros([self._Np,3])
        i = 0
        while i < self._Np:            
            coord = np.array([np.random.uniform(0,self._Lx,1),np.random.uniform(0,self._Ly,1),np.random.uniform(0,self._Lz,1)]).T  
            if self._reject(coord) == False:
                coords[i]=coord
                i += 1 
        #Seeding Code
        #Uniform Random Generator
        #coords = np.array([np.random.uniform(0,self._Lx,self._Np),np.random.uniform(0,self._Ly,self._Np),np.random.uniform(0,self._Lz,self._Np)]).T
        
        self['pore.coords'] = coords
        self._logger.debug(sys._getframe().f_code.co_name+": End of method") 
        
    def _generate_cubic_pores(self):
        self._logger.info(sys._getframe().f_code.co_name+": Place pores on cubic grid and jiggle")
        Nx = 10
        Ny = 10
        Nz = 10
        self._Np = Nx*Ny*Nz
        equal_spacing = cbrt(self._Np)*3e-6
        ind = np.arange(0,self._Np)
        aniso_spacing = equal_spacing*self._scale_factor
        Px = aniso_spacing[0]*(0.5+np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=np.float).T)[:,0]+aniso_spacing[0]*(np.random.rand(self._Np)-0.5)/5
        Py = aniso_spacing[1]*(0.5+np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=np.float).T)[:,1]+aniso_spacing[0]*(np.random.rand(self._Np)-0.5)/5
        Pz = aniso_spacing[2]*(0.5+np.array(np.unravel_index(ind, dims=(Nx, Ny, Nz), order='F'),dtype=np.float).T)[:,2]+aniso_spacing[0]*(np.random.rand(self._Np)-0.5)/5
        coords = np.vstack((Px,Py,Pz)).T
        self.set_pore_data(prop='coords',data=coords)
        self.set_pore_info(label='all',locations=np.ones_like(coords[:,0]))
        self._logger.debug(sys._getframe().f_code.co_name+": End of method")
        self._Lx = Px.max()+aniso_spacing[0]/2
        self._Ly = Py.max()+aniso_spacing[1]/2
        self._Lz = Pz.max()+aniso_spacing[2]/2
    
    def _prob_func(self,m):
        a = 35
        b = 0.2
        p = ((m**a) + ((1-m)**a) + (2*b))/(1 + (2*b))
        #p = ((m**a) + b)/(1 + b)
        return p

    def _reject(self,point):
    
        x = point[0,0]
        y = point[0,1]
        #z = point[0,2]
        Px = self._prob_func(x)
        Py = self._prob_func(y)
        #Pz = prob_func(z)
        nrand = np.random.uniform(0,1,1)
    
        if Px < nrand and Py < nrand:
            rejection = True
        else:
            rejection = False
    
        return rejection
    
    def _generate_throats(self):
        r"""
        Generate the throats connections
        """
        self._logger.info(sys._getframe().f_code.co_name+": Define connections between pores")
        Np = self._Np
        pts = self['pore.coords']
        #Generate 6 dummy domains to pad onto each face of real domain
        #This prevents surface pores from making long range connections to each other
        Lx = self._Lx
        Ly = self._Ly
        Lz = self._Lz

        #Reflect in X = Lx and 0
        Pxp = pts.copy()
        Pxp[:,0]=(2*Lx-Pxp[:,0])
        Pxm= pts.copy()
        Pxm[:,0] = Pxm[:,0]*(-1)
        #Reflect in Y = Ly and 0
        Pyp = pts.copy()
        Pyp[:,1]=(2*Ly-Pxp[:,1])
        Pym = pts.copy()
        Pym[:,1] = Pxm[:,1]*(-1)
        #Reflect in Z = Lz and 0
        Pzp = pts.copy()
        Pzp[:,2]=(2*Lz-Pxp[:,2])
        Pzm = pts.copy()
        Pzm[:,2] = Pxm[:,2]*(-1)
        #Add dummy domains to real domain
        pts = np.vstack((pts,Pxp,Pxm,Pyp,Pym,Pzp,Pzm)) #Order important for boundary logic
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
        self['throat.conns']=sp.vstack((adjmat.row, adjmat.col)).T
        self['pore.all'] = np.ones(len(self['pore.coords']), dtype=bool)
        self['throat.all'] = np.ones(len(self['throat.conns']), dtype=bool)

        #New code to identify boundary pores - those that connect to pores inside and outside original set of pores 
        boundary_pore_list = []
        for i in sp.arange(0,sp.shape(Tri.simplices)[0]):
            pores_in = Tri.simplices[i] < Np # Pores in the original domain
            if (sum(pores_in) >= 1) and (sum(pores_in) < len(pores_in)):
                for j in range(len(Tri.simplices[i])):
                    if pores_in[j] == True:
                        pore_id = Tri.simplices[i][j]
                        if pore_id not in boundary_pore_list:
                            boundary_pore_list.append(pore_id)
        
        vor_bounds=sp.asarray(boundary_pore_list)

        self.set_info(pores=vor_bounds,label='inner_boundary')
        self.set_info(pores='all',label='internal')
        self.set_info(throats='all',label='internal')
        # Do Voronoi diagram - creating voronoi polyhedra around each pore and save vertex information
        vor = Voronoi(pts)
        all_verts = sp.ndarray(Np,dtype=object)
        for i,polygon in enumerate(vor.point_region[0:Np]):
            if -1 not in vor.regions[polygon]:
                all_verts[i]=np.around(vor.vertices[vor.regions[polygon]],10)
            else:
                all_verts[i]="unbounded"
        self['pore.vertices']=all_verts
        self._logger.debug(sys._getframe().f_code.co_name+": End of method")
        #Add new pores at external throat centers to create coplanar boundaries
        self.boundary_pores()
        self["pore.coords"]=self["pore.coords"]*self._rescale_factor
        for i in range(len(self["pore.vertices"])):
            self["pore.vertices"][i]=self["pore.vertices"][i]*self._rescale_factor
        external_pores = self.pores(labels='internal',mode='difference')
        external_throats = self.throats(labels='internal',mode='difference')
        self.set_info(pores=external_pores,label='external')
        self.set_info(throats=external_throats,label='external')


    def _add_labels(self):
        r'''
        This finds surface pores simply by proximity to the domain boundaries.
        A better approach is necessary 
        '''
        coords = self.get_pore_data(prop='coords')
        self.set_pore_info(label='front',locations=(coords[:,0]<(0.1*self._Lx)))
        self.set_pore_info(label='back',locations=(coords[:,0]>(0.9*self._Lx)))
        self.set_pore_info(label='left',locations=(coords[:,1]<(0.1*self._Ly)))
        self.set_pore_info(label='right',locations=(coords[:,1]>(0.9*self._Ly)))
        self.set_pore_info(label='bottom',locations=(coords[:,2]<(0.1*self._Lz)))
        self.set_pore_info(label='top',locations=(coords[:,2]>(0.9*self._Lz)))
        bnds = self.pores(labels=['front','back','left','right','bottom','top'])
        self.set_pore_info(label='boundary',locations=bnds)
        
        
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
        pnum = self.pores()
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
    
    def boundary_pores(self):
        r"""
        This method runs through the boundary pores and identifies the throats 
        that align with the boundary plane. As there are no connections to the 
        dummy pores we also need to generate the pores at the centroids of the 
        boundary throats or at the boundary coordinate associated with the 
        plane leaving the other pore coords the same.  This will mean that 
        normal vectors are correct for throat props.
        """
        boundary_pores = self.pores('inner_boundary')
        '''Traditionally x will be back and front, y is left and right and z is 
        top and bottom.  However this could change in future and labels might 
        also change so best to do things from scratch.'''
        #Look at the number of occurences of a coordinate for the pore's vertices, if 3 or more we have a face
        boundary_throats = []
        new_boundary_pores = []
        throat_centers = []
        new_conns=[]
        pore_coords=[]
        #Find boundary extent
        [x_min,x_max,y_min,y_max,z_min,z_max]=self.vertex_dimension(self.pores(),parm='minmax')
        min_point = np.around(np.array([x_min,y_min,z_min]),10)
        max_point = np.around(np.array([x_max,y_max,z_max]),10)
        delta = (max_point - min_point)*0.25
        N = self.num_pores()
        new_throat_count = 0
        for pore in boundary_pores:
            verts = np.around(self["pore.vertices"][pore],10)
            #Cycle through coordinates
            for i in range(3):
                throat_verts = []
                extruded_verts = []
                temp_coord = sp.copy(self["pore.coords"][pore])
                freq = st.itemfreq(verts[:,i])
                for coord in freq:
                    #if more than 2 occurences
                    if coord[1]>2:
                        #Pick up planar value for referencing
                        value = coord[0]
                        #Are we at minimum or maximum extent?
                        #Extrude into the domain so as not to alter total dimensions
                        extrude_value = np.zeros(3)
                        if value == min_point[i]:
                            extrude_value[i] = delta[i]
                            temp_coord[i] = min_point[i]
                        else:
                            extrude_value[i] = -delta[i]
                            temp_coord[i] = max_point[i]
                        for vert in verts:
                            #Pick up all verts in the plane to define the throat
                            if vert[i] == value:
                                throat_verts.append(vert)
                                extruded_verts.append(vert+extrude_value)
                                
                #If we found a planar throat then add to the list
                if len(throat_verts) > 0:
                    new_conns.append(np.array([pore,new_throat_count+N]))
                    new_throat_count += 1
                    new_boundary_pores.append(np.asarray(throat_verts+extruded_verts))
                    throat_verts = np.asarray(throat_verts)
                    boundary_throats.append(throat_verts)
                    throat_centers.append(sp.array([throat_verts[:,0].mean(),throat_verts[:,1].mean(),throat_verts[:,2].mean()]))
                    pore_coords.append(temp_coord)
            
        #First Attempt to try and create new pores and new connections
        #Add new pores and connections
        self.extend(pore_coords=pore_coords, throat_conns=new_conns)
        #Record new number of pores
        M = self.num_pores()
        new_pores = np.arange(N,M)
        #Identify which boundary the pore sits on
        front = self.pores()[self['pore.coords'][:,0]==min_point[0]]
        back = self.pores()[self['pore.coords'][:,0]==max_point[0]]
        left = self.pores()[self['pore.coords'][:,1]==min_point[1]]
        right = self.pores()[self['pore.coords'][:,1]==max_point[1]]
        bottom = self.pores()[self['pore.coords'][:,2]==min_point[2]]
        top = self.pores()[self['pore.coords'][:,2]==max_point[2]]
        #Assign labels
        self.set_info(pores=new_pores,label='boundary')        
        self.set_info(pores=right,label='right_boundary') 
        self.set_info(pores=left,label='left_boundary') 
        self.set_info(pores=front,label='front_boundary') 
        self.set_info(pores=back,label='back_boundary') 
        self.set_info(pores=top,label='top_boundary') 
        self.set_info(pores=bottom,label='bottom_boundary')
        #Apply the extruded throat verts and original boundary throat verts to create new pore volume
        #These volumes may need some attention later
        self["pore.vertices"][N:M] = new_boundary_pores
		
    def vertex_dimension(self,face1=[],face2=[],parm='volume'):
        r"""
        Return the domain extent based on the vertices
        
        This function is better than using the pore coords as they may be far 
        away from the original domain size.  And will alter the effective 
        properties which should be based on the original domain sizes. Takes 
        one or two sets of pores and works out different geometric properties
        if "length" is specified and two lists are given the planarity is 
        determined and the appropriate length (x,y,z) is returned.  It should 
        work the same as domain length and area if vertices are not in network 
        by using coordinates.
        
        e.g.    vertex_extent(face1=inlet,face2=outlet,parm='volume')
                vertex_extent(geom.pores(),parm='area_xy')
                vertex_extent(face1=inlet,parm='area')
                vertex_extent(face1=inlet,face2=outlet,parm='length')

        """
        pores=np.array([],dtype=int)
        if len(face1)>0:
            pores=np.hstack((pores,face1))
        if len(face2)>0:
            pores=np.hstack((pores,face2))
        
        face1_coords = self["pore.coords"][face1]
        face2_coords = self["pore.coords"][face2]
        face1_planar = np.zeros(3)
        face2_planar = np.zeros(3)
        planar = np.zeros(3)
        for i in range(3):
            if len(np.unique(face1_coords[:,i]))==1:
                face1_planar[i]=1
            if len(np.unique(face2_coords[:,i]))==1:
                face2_planar[i]=1
        if len(face1)>0 and len(face2)>0:
            planar = face1_planar*face2_planar
        elif len(face1)>0:
            planar = face1_planar
        elif len(face2)>0:
            planar = face2_planar
        else:
            return 0
            
        if "pore.vertices" in self.props():
            vert_list = self["pore.vertices"][pores]
        else:
            vert_list = self["pore.coords"][pores]
        vx_min = 1e32
        vx_max = -1e32
        vy_min = 1e32
        vy_max = -1e32
        vz_min = 1e32
        vz_max = -1e32
        output = 0
        for verts in vert_list:
            if verts[:,0].min()<vx_min:
                vx_min=verts[:,0].min()
            if verts[:,0].max()>vx_max:
                vx_max=verts[:,0].max()
            if verts[:,1].min()<vy_min:
                vy_min=verts[:,1].min()
            if verts[:,1].max()>vy_max:
                vy_max=verts[:,1].max()
            if verts[:,2].min()<vz_min:
                vz_min=verts[:,2].min()
            if verts[:,2].max()>vz_max:
                vz_max=verts[:,2].max()
        width = np.around(vx_max-vx_min,10)
        depth = np.around(vy_max-vy_min,10)
        height = np.around(vz_max-vz_min,10)
        if parm == 'volume':
            output =  width*depth*height
        elif parm == 'area_xy' or (parm == 'area' and planar[2]==1):
            output = width*depth
        elif parm == 'area_xz' or (parm == 'area' and planar[1]==1):
            output = width*height
        elif parm == 'area_yz' or (parm == 'area' and planar[0]==1):
            output = depth*height
        elif parm == 'length_x' or (parm == 'length' and planar[0]==1):
            output = width
        elif parm == 'length_y'or (parm == 'length' and planar[1]==1):
            output = depth
        elif parm == 'length_z'or (parm == 'length' and planar[2]==1):
            output = height
        elif parm == 'minmax':
            output = [vx_min,vx_max,vy_min,vy_max,vz_min,vz_max]
        return output
    
    def domain_length(self,face_1,face_2):
        r"""
        Returns the distance between two faces
        No coplanar checking this is done in vertex_dimension
        """
        L = self.vertex_dimension(face_1,face_2,parm='length')
        return L

    def domain_area(self,face):
        r"""
        Returns the area of a face
        No coplanar checking this is done in vertex_dimension
        """
        A = self.vertex_dimension(face,parm='area')
		
        return A       

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
