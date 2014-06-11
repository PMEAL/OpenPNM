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
import _transformations as tr
from scipy.spatial import ConvexHull
from math import atan2

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
    >>> pn.generate(num_pores=100,domain_size=[100,100,100])
    >>> pn.num_pores()
    100
    >>> type(pn.num_throats())
    <class 'numpy.int32'>

    """

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
        self._generate_pores()
        self._generate_throats()
#        self._add_boundaries()
        self._add_labels()
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

    def _centroid(self,points):
        import numpy as np
        columns = np.shape(points)[1]
        if columns == 2:
            cen = np.array([points[:,0].mean(),points[:,1].mean()])
        elif columns == 3:
            cen = np.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
        else:
            cen = 0
        return cen    
    
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
        pts = np.vstack((pts,Pxp,Pxm,Pyp,Pym,Pzp,Pzm))
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
        
        # Do Voronoi diagram - creating voronoi polyhedra around each pore and save vertex information
        vor = Voronoi(pts)
        all_verts = sp.ndarray(Np,dtype=object)
        centroids = sp.ndarray(Np,dtype=object)
        for i,polygon in enumerate(vor.point_region[0:Np]):
            if -1 not in vor.regions[polygon]:
                all_verts[i]=vor.vertices[vor.regions[polygon]]
                centroids[i]=self._centroid(all_verts[i])
            else:
                all_verts[i]="unbounded"
                centroids[i]=pts[i]
        self.set_pore_data(prop='vertices',data=all_verts)
        self.set_pore_data(prop='centroid',data=centroids)
        self._add_throat_props()
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
    
    def _add_throat_props(self,radius=1e-06):
        r"""
        This method does all the throat properties for the voronoi cages 
        including offseting the vertices surrounding each pore by an amount that
        replicates erroding the facet of each throat by the fibre radius 

        For each connection or throat find the shared vertices
        Rotate the vertices to align with the xy-plane and get rid of z-coordinate
        Merge vertices within a certain distance of each other based on the range of data points
        Compute the convex hull of the 2D points giving a set of simplices which define neighbouring vertices in a clockwise fashion
        For each triplet calculate the offset position given the fibre radius
        Translate back into 3D
        """
        connections = self.get_throat_data(prop='conns')
        coords = self.get_pore_data(prop='coords')
        verts = self.get_pore_data(prop='vertices')
        normals = coords[connections[:,0]]-coords[connections[:,1]]
        area = sp.ndarray(len(connections),dtype=object)
        perimeter = sp.ndarray(len(connections),dtype=object)
        centroids = sp.ndarray(len(connections),dtype=object)
        offset_verts = sp.ndarray(len(connections),dtype=object)
        for i,throat_pair in enumerate(connections):
            shared_verts = []
            pore_a = throat_pair[0]
            pore_b = throat_pair[1]
            " Identify shared verts "
            for vert_a in verts[pore_a]:
                for vert_b in verts[pore_b]:
                    if (vert_a[0] == vert_b[0]) and (vert_a[1] == vert_b[1]) and (vert_a[2] == vert_b[2]):
                        shared_verts.append(vert_a)
            if len(shared_verts) >=3:
                shared_verts = np.asarray(shared_verts)
                area[i],perimeter[i],offset_verts[i] = self._get_throat_geom(shared_verts,normals[i])
                centroids[i]=self._centroid(shared_verts)
            else:
                area[i]=0.0

        self.set_data(prop='area',throats='all',data=area)
        self.set_data(prop='perimeter',throats='all',data=perimeter)
        self.set_throat_data(prop='centroid',data=centroids)
        self.set_throat_data(prop='offset_verts',data=offset_verts)
        " Temporary Code to remove the smallest 2% of throat area connections "
        " This can be used in algorithms to ignore certain connections if required - like in the range of capillary pressures in OP "
        #smallest_throat = min(area)
        #largest_throat = max(area)
        average_area = sp.mean(area)
        cutoff = average_area/100
        #remove_range = smallest_throat + ((largest_throat-smallest_throat)*0.02)
        excluded_throats = []
        for i,throat in enumerate(connections):
            if area[i]<=cutoff:
            #if area[i]==0.0:
                excluded_throats.append(i)
        excluded_throats = np.asarray(excluded_throats)
        if len(excluded_throats) > 0:
            self.trim(throats=excluded_throats)
    
    def _get_throat_geom(self,verts,normal):
        z_axis = [0,0,1]
        " For boundaries some facets will already be aligned with the axis - if this is the case a rotation is unnecessary and could also cause problems "
        angle = tr.angle_between_vectors(normal,z_axis)
        if (angle==0.0)or(angle==np.pi):
            "We are already aligned"
            rotate_input = False
            facet = verts
        else:
            rotate_input = True
            M = tr.rotation_matrix(tr.angle_between_vectors(normal,z_axis),tr.vector_product(normal,z_axis))
            facet = np.dot(verts,M[:3,:3].T)
        x = facet[:,0]
        y = facet[:,1]
        z = facet[:,2]
        if (np.around(z.std(),3)!=0.000):
            print("Rotation failed")
        facet_coords_2D = np.column_stack((x,y))
        hull = ConvexHull(facet_coords_2D)
        verts_2D = facet_coords_2D[hull.vertices]
        fused_verts=self.fuse_verts(verts_2D)
        if len(fused_verts) <3:
            #we fused too many
            fused_verts=self.fuse_verts(verts_2D,0.0025)
        offset = []
        fibre_rad = 3e-06
        for i,vert in enumerate(fused_verts):
            " Collect three adjacent points and compute the offset of the first "
            triplet = (vert, np.roll(fused_verts,-1,axis=0)[i],np.roll(fused_verts,1,axis=0)[i])
            offset.append(self.offset_vertex(triplet,fibre_rad))
        offset = np.asarray(offset)
    
        original_area = self.PolyArea2D(verts_2D)            
        all_points = np.concatenate((verts_2D,offset),axis=0)
        try:
            total_hull = ConvexHull(all_points,qhull_options='Pp') #ignores very small angles
            total_area = self.PolyArea2D(all_points[total_hull.vertices])
        except sp.spatial.qhull.QhullError:
            print(all_points)
            total_area =999

        if (total_area>original_area): # Throat is fully occluded
            throat_area = 0.0
            throat_perimeter = 0.0
            output_offset = []
        else:
            offset_hull = ConvexHull(offset)
            offset_verts_2D = offset[offset_hull.vertices]
            throat_area = self.PolyArea2D(offset_verts_2D)
            throat_perimeter = self.PolyPerimeter2D(offset_verts_2D)
            " Make 3D again in rotated plane "
            offset_verts_3D = np.column_stack((offset_verts_2D,z[0:len(offset_verts_2D)]))
            " Get matrix to un-rotate the co-ordinates back to the original orientation if we rotated in the first place"
            if (rotate_input):
                M1 = tr.inverse_matrix(M)
                " Unrotate the offset coordinates "
                output_offset = np.dot(offset_verts_3D,M1[:3,:3].T)
            else:
                output_offset = offset_verts_3D

        return throat_area, throat_perimeter, output_offset
    
    def offset_vertex(self,points,rad = 0.01):
        " We are passed in a set of 3 points forming vertices of two adjoining simplexes of the convex hull of a voronoi facet "
        " We need to offset the vertices normal to the fibre direction (or adjoining vectors) by the fibre radius "
        " This is achieved by finding the half angle between the two adjoining vectors and a direction "
        " Mid-point must be the first in the array "
        p0 = np.array(points[0])
        p1 = np.array(points[1])
        p2 = np.array(points[2])
        " Now make the midpoint the origin "
        vector1 = p1-p0
        vector2 = p2-p0

        " Find what quadrant the vector is pointing in - atan2 function takes account of signs "
        " 0 means aligned with x-axis, pi is aligned with -xaxis, positive numbers are positive y and negative numbers are negative y "
        " The angle between the vectors should always be within 180 degrees of one another in a convex hull "

        q1 = atan2(vector1[1],vector1[0])
        q2 = atan2(vector2[1],vector2[0])
        alpha = 0.5*tr.angle_between_vectors(vector1,vector2)
    
        " We always want to offset from the first vertex we get to - going anti-clockwise from the x-axis "
        " Check if both vectors point up or both point down - if so the first one we get to will have smaller q value "
        if q1*q2 >=0.0:
            if q1<q2:
                theta = q1
            else:
                theta = q2
        else:
            "if vector 1 is more rotated away from positive xaxis than vector 2 is rotated away from negative xaxis - use it"
            " and vice-versa "
            if (abs(q1)+abs(q2)>np.pi):
                "vectors are pointing negative x so take whichever has positive q-value - like a pacman facing left"
                if q1>=0:
                    theta = q1
                else:
                    theta = q2
            else:
                "vectors are pointing positive x so take whichever is negative"
                if q1<=0:
                    theta = q1
                else:
                    theta = q2
    
        x = rad*np.cos(alpha+theta)/np.sin(alpha)
        y = rad*np.sin(alpha+theta)/np.sin(alpha)

        "Add the midpoint back in"
        output = [x+p0[0],y+p0[1]]    
    
        return output

    def dist2(self,p1, p2):
        return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2

    def fuse(self,points, d):
        ret = []
        d2 = d * d
        n = len(points)
        taken = [False] * n
        for i in range(n):
            if not taken[i]:
                count = 1
                point = [points[i][0], points[i][1]]
                taken[i] = True
                for j in range(i+1, n):
                    if self.dist2(points[i], points[j]) < d2:
                        point[0] += points[j][0]
                        point[1] += points[j][1]
                        count+=1
                        taken[j] = True
                point[0] /= count
                point[1] /= count
                ret.append((point[0], point[1]))
        return ret

    def fuse_verts(self,verts,percentage=0.05):
        #Work out largest span
        x_span = max(verts[:,0])- min(verts[:,0])
        y_span = max(verts[:,1])- min(verts[:,1])
        if x_span > y_span:
            tolerance = x_span*percentage
        else:
            tolerance = y_span*percentage
        #fuse vertices lying within 5% of the largest span    
        return self.fuse(verts,tolerance)

    def PolyArea2D(self,pts):
        lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
        area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
        return area

    def PolyPerimeter2D(self,pts):
        lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
        perimeter = sum(np.sqrt((x2-x1)**2+(y2-y1)**2) for x1,y1,x2,y2 in lines)
        return perimeter

if __name__ == '__main__':
    #Run doc tests
    import doctest
    doctest.testmod(verbose=True)
