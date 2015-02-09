r"""
===============================================================================
pore_volume --
===============================================================================

"""
import scipy as _sp
import numpy as np
from scipy.spatial import Delaunay
from scipy.spatial import ConvexHull
import OpenPNM.Utilities.misc as misc
from multiprocessing import Pool
from scipy import ndimage
from scipy.io import savemat
import gc

def _get_vertex_range(verts):
    "Find the extent of the vetrices"
    [vxmin,vxmax,vymin,vymax,vzmin,vzmax]=[1e20,0,1e20,0,1e20,0]
    for vert in verts:
        if np.min(vert[:,0])< vxmin:
            vxmin = np.min(vert[:,0])
        if np.max(vert[:,0])> vxmax:
            vxmax = np.max(vert[:,0])
        if np.min(vert[:,1])< vymin:
            vymin = np.min(vert[:,1])
        if np.max(vert[:,1])> vymax:
            vymax = np.max(vert[:,1])
        if np.min(vert[:,2])< vzmin:
            vzmin = np.min(vert[:,2])
        if np.max(vert[:,2])> vzmax:
            vzmax = np.max(vert[:,2])
    return [vxmin,vxmax,vymin,vymax,vzmin,vzmax]
    
def _get_fibre_image(network,cpores,vox_len,fibre_rad,verts,export_mat=False,mat_file='OpenPNMFibres'):
    r"""
    Produce image by filling in voxels along throat edges using Bresenham line
    Then performing distance transform on fibre voxels to erode the pore space
    """

    cthroats = network.find_neighbor_throats(pores=cpores)
    #geom_throats = network.map_throats(geometry,cthroats,return_mapping=True)["target"]
    
    cverts = verts[cthroats]
    [vxmin,vxmax,vymin,vymax,vzmin,vzmax] = _get_vertex_range(cverts)
    "Translate vertices so that minimum occurs at the origin"
    for index in range(len(cverts)):
        cverts[index] -= np.array([vxmin,vymin,vzmin])
    "Find new size of image array"
    cdomain=np.around(np.array([(vxmax-vxmin),(vymax-vymin),(vzmax-vzmin)]),6)
    print("Fibre Domain Range: "+str(np.around(cdomain,5)))            
    lx = np.int(np.around(cdomain[0]/vox_len)+1)
    ly = np.int(np.around(cdomain[1]/vox_len)+1)
    lz = np.int(np.around(cdomain[2]/vox_len)+1)
    "Define voxel image of pore space"
    pore_space=np.ndarray([lx,ly,lz],dtype=np.uint8)
    pore_space.fill(1)
    "Get image of the fibres"
    line_points = bresenham(cverts,vox_len)
    line_ints = (np.around(line_points/vox_len,0)).astype(int)
    for x,y,z in line_ints:
        pore_space[x][y][z]=0
    #print("about to dt")
    chunk_len = 100
    if (lx > chunk_len) or (ly > chunk_len) or (lz > chunk_len):
        cx = np.ceil(lx/chunk_len).astype(int)
        cy = np.ceil(ly/chunk_len).astype(int)
        cz = np.ceil(lz/chunk_len).astype(int)
    else:
        cx = cy = cz = 1
    fibre_space = np.ndarray(shape=[lx,ly,lz],dtype=np.uint8)
    for ci in range(cx):
        for cj in range(cy):
            for ck in range(cz):
                "Work out chunk range"
                cxmin = ci*chunk_len
                cxmax = (ci+1)*chunk_len + 2*fibre_rad.astype(int)
                cymin = cj*chunk_len
                cymax = (cj+1)*chunk_len + 2*fibre_rad.astype(int)
                czmin = ck*chunk_len
                czmax = (ck+1)*chunk_len + 2*fibre_rad.astype(int)
                "Don't overshoot"
                if cxmax > lx:
                    cxmax = lx
                if cymax > ly:
                    cymax = ly
                if czmax > lz:
                    czmax = lz
                dt = ndimage.distance_transform_edt(pore_space[cxmin:cxmax,cymin:cymax,czmin:czmax])
                fibre_space[cxmin:cxmax,cymin:cymax,czmin:czmax][dt<=fibre_rad]=0
                fibre_space[cxmin:cxmax,cymin:cymax,czmin:czmax][dt>fibre_rad]=1
                del dt
    del pore_space
    if export_mat:
        matlab_dict = {"fibres":fibre_space}
        savemat(mat_file+str(lx)+str(ly)+str(lz),matlab_dict,format='5',long_field_names=True)
    porosity = np.around(np.sum(fibre_space)/np.size(fibre_space),3)
    print("Porosity from fibre image: "+str(porosity) )
    print("Size of Fibre Space: "+str(np.size(fibre_space)))
    return fibre_space

def _get_voxel_volume(network,chunk,vox_len,fibre_rad,verts,fibre_image):
    r"""
    Calculate the volume by divinding space into voxels, working out nearest neighours to get hulls
    returns number of voxels in pore both fibre and open space portions
    """
    "Unpack the pores and range of the chunk"
    [cpores,crange]=chunk
    [cxmin,cxmax,cymin,cymax,czmin,czmax]=crange
    #network,cpores,vox_len,fibre_rad,verts = chunk_data
    voxel = vox_len**3
    
    "Boundary pores have zero volume so we ignore them, only evaluate nearest neighbors"
    "For points directly in chunk and those connected by network"
    nbps = network.pores('boundary',mode='not')#get list of non-boundary pores, these have zero volume
    neighbors = network.find_neighbor_pores(pores=cpores)
    my_pores = np.asarray(list(set(np.concatenate((cpores,neighbors))).intersection(set(nbps))))
    if len(my_pores)==0:
        return np.zeros(len(cpores))
    my_points = network["pore.coords"][my_pores]
    pore_volume = _sp.zeros(len(my_pores),dtype=int)
    fibre_volume = _sp.zeros(len(my_pores),dtype=int) # needed for when domain is compressed or stretched to ensure total fibre vol conserved

    stp = np.array([cxmin,cymin,czmin]) #starting point
    "Translate points so starting point of chunk is at origin"
    temp_points = my_points - stp
    si = np.around(stp/vox_len).astype(int) # starting index
    "Find new size of image array"
    cdomain=np.array([(cxmax-cxmin),(cymax-cymin),(czmax-czmin)])
    print("Vertex Range: "+str(np.around(crange,5)))            
    lx = np.int(np.around(cdomain[0]/vox_len)+1)
    ly = np.int(np.around(cdomain[1]/vox_len)+1)
    lz = np.int(np.around(cdomain[2]/vox_len)+1)
    "Get chunk of fibre image"
    fibre_space = fibre_image[si[0]:si[0]+lx,si[1]:si[1]+ly,si[2]:si[2]+lz]
    
    "Hull Method 5 - Skilearn Neighbors"
    hull_space=np.zeros([lx,ly,lz],dtype=np.uint16)
    from sklearn.neighbors import NearestNeighbors
    temp_points /= vox_len
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(temp_points)
    qs = []
    for i in range(lx):
        for j in range(ly):
            for k in range(lz):
                qs.append([i,j,k])
    indices = nbrs.kneighbors(qs,return_distance=False)
    for n, [i,j,k] in enumerate(qs):
        hull_space[i][j][k]=my_pores[indices[n]]
    del indices
    del nbrs
            
    for index,pore in enumerate(my_pores):
        pore_volume[index] = np.sum((fibre_space == 1)&(hull_space==pore))
        fibre_volume[index] = np.sum((fibre_space == 0)&(hull_space==pore))
    print("Size of Chunk Space: "+str(np.size(hull_space)))
    del fibre_space
    del hull_space
    gc.collect()
    
    return pore_volume,fibre_volume,my_pores
    
def bresenham(faces,dx):
    line_points=[]
    for face in faces:
        #Get in hull order
        f2d = face[:,0:2]
        hull = ConvexHull(f2d,qhull_options='QJ Pp')
        face=face[hull.vertices]
        for i in range(len(face)):
            vec = face[i]-face[i-1]
            vec_length = np.linalg.norm(vec)
            increments = np.ceil(vec_length/dx)
            check_p_old = np.array([-1,-1,-1])
            for x in np.linspace(0,1,increments):
                check_p_new = face[i-1]+(vec*x)
                if np.sum(check_p_new - check_p_old) != 0:
                    line_points.append(check_p_new)
                    check_p_old = check_p_new
    return np.asarray(line_points)

def _get_hull_volume(points):
    r"""
    Calculate the volume of a set of points by dividing the bounding surface into triangles and working out the volume of all the pyramid elements
    connected to the volume centroid
    """
    " remove any duplicate points - this messes up the triangulation "
    points = _sp.asarray(misc.unique_list(np.around(points,10)))
    try:
        tri = Delaunay(points,qhull_options='QJ Pp')
    except _sp.spatial.qhull.QhullError:
        print(points)
    " We only want points included in the convex hull to calculate the centroid "
    hull_centroid = _sp.array([points[:,0].mean(),points[:,1].mean(),points[:,2].mean()])
    hull_volume = 0.0
    pyramid_COMs = []
    for ia, ib, ic in tri.convex_hull:
        " Points making each triangular face "
        " Collection of co-ordinates of each point in this face "
        face_x = points[[ia,ib,ic]][:,0]
        face_y = points[[ia,ib,ic]][:,1]
        face_z = points[[ia,ib,ic]][:,2]
        " Average of each co-ordinate is the centroid of the face "
        face_centroid = [face_x.mean(),face_y.mean(),face_z.mean()]
        face_centroid_vector = face_centroid - hull_centroid
        " Vectors of the sides of the face used to find normal vector and area "
        vab = points[ib] - points[ia]
        vac = points[ic] - points[ia]
        vbc = points[ic] - points[ib] # used later for area
        #face_COM = (vab+vac)/3
        " As vectors are co-planar the cross-product will produce the normal vector of the face "
        face_normal = _sp.cross(vab,vac)
        try:
            face_unit_normal = face_normal/_sp.linalg.norm(face_normal)
        except RuntimeWarning:
            print("Pore Volume Error:" +str(vab)+" "+str(vac))
        " As triangles are orientated randomly in 3D we could either transform co-ordinates to align with a plane and perform 2D operations "
        " to work out the area or we could work out the lengths of each side and use Heron's formula which is easier"
        " Using Delaunay traingulation will always produce triangular faces but if dealing with other polygons co-ordinate transfer may be necessary "
        a = _sp.linalg.norm(vab)
        b = _sp.linalg.norm(vbc)
        c = _sp.linalg.norm(vac)
        " Semiperimeter "
        s = 0.5*(a+b+c)
        face_area = _sp.sqrt(s*(s-a)*(s-b)*(s-c))
        " Now the volume of the pyramid section defined by the 3 face points and the hull centroid can be calculated "
        pyramid_volume = _sp.absolute(_sp.dot(face_centroid_vector,face_unit_normal)*face_area/3)
        " Each pyramid is summed together to calculate the total volume "
        hull_volume += pyramid_volume
        " The Centre of Mass will not be the same as the geometrical centroid "
        " A weighted adjustment can be calculated from the pyramid centroid and volume "
        vha = points[ia]-hull_centroid
        vhb = points[ib]-hull_centroid
        vhc = points[ic]-hull_centroid
        pCOM = ((vha+vhb+vhc)/4)*pyramid_volume
        pyramid_COMs.append(pCOM)
    if _sp.isnan(hull_volume):
        hull_volume = 0.0
    if hull_volume>0:
        hull_COM = hull_centroid + _sp.mean(_sp.asarray(pyramid_COMs),axis=0)/hull_volume
    else:
        hull_COM = hull_centroid

    return hull_volume, hull_COM

def sphere(geometry,
           pore_diameter='pore.diameter',
           **kwargs):
    r"""
    Calculate pore volume from diameter for a spherical pore body
    """
    diams = geometry[pore_diameter]
    value=_sp.pi/6*diams**3
    return value

def cube(geometry,
         pore_diameter='pore.diameter',
         **kwargs):
    r"""
    Calculate pore volume from diameter for a cubic pore body
    """
    diams = geometry[pore_diameter]
    value = diams**3
    return value

def voronoi(network,
            geometry,
            **kwargs):
    r"""
    Calculate volume from the convex hull of the offset vertices making the throats surrounding the pore
    Also calculate the centre of mass for the volume
    """
    pores = geometry.map_pores(network,geometry.pores())
    Np = len(pores)
    volume = _sp.zeros(Np)
    com = _sp.zeros([Np,3])
    for i in range(Np):
        throat_vert_list = []
        net_throats=network.find_neighbor_throats([pores[i]])
        geom_throats = network.map_throats(target=geometry,throats=net_throats,return_mapping=True)['target']
        if len(geom_throats) > 1:
            for throat in geom_throats:
                geom_throat_verts = geometry["throat.offset_vertices"][throat]
                if geom_throat_verts is not None:
                    for j in range(len(geom_throat_verts)):
                        throat_vert_list.append(geom_throat_verts[j])
            throat_array=_sp.asarray(throat_vert_list)
            if len(throat_array)>4:
                volume[i],com[i] = _get_hull_volume(throat_array)
            else:
                volume[i]=0
        elif len(geom_throats) == 1 and 'throat.centroid' in geometry.props():
                com[i]=geometry['throat.centroid'][geom_throats]
                volume[i]=0
    "Find any pores with centroids at origin and use the mean of the pore vertices instead"
    "Not doing this messes up hydraulic conductances using centre to centre"
    ps = np.where(~com.any(axis=1))[0]
    if len(ps) >0:
        for pore in ps:
            com[pore]=np.mean(geometry["pore.vertices"][pore],axis=0)
    geometry["pore.centroid"]=com

    return volume

def voronoi_vox(network,
                geometry,
                fibre_rad,
                export_mat='False',
                mat_file='mat_file',
                **kwargs):
    r"""
    Compute the pore volumes by creating a voxel image of the domain with Bresenham lines creating fibres
    Warning this method takes longer than standard voronoi calculation
    Maximum domain size for one chunk is 700**3
    N.B may be inefficient as whole network is calculated and then specific geometry returned
    """
    
    import OpenPNM.Utilities.vertexops as vo
    #from scipy.io import savemat
    
    Np = network.num_pores()
    geom_pores = geometry.map_pores(network,geometry.pores())
    volume = _sp.zeros(Np)
    fibre_vox = _sp.zeros(Np,dtype=np.uint16)
    vox_len=1e-6
    voxel = vox_len**3

    B1 = network.pores("left_boundary")
    B2 = network.pores("right_boundary")
    [xmin,xmax,ymin,ymax,zmin,zmax]=vo.vertex_dimension(network,B1,B2,parm='minmax')
    domain=np.array([(xmax-xmin),(ymax-ymin),(zmax-zmin)])    
    fibre_rad = np.around((fibre_rad-(vox_len/2))/vox_len,0).astype(int) #voxel length

    points = network["pore.coords"]
    #verts = network["throat.vertices"]
    "Below method copied from geometry model throat.vertices"
    "Needed now as network may not have all throats assigned to geometry"
    "i.e network['throat.vertices'] could return garbage"
    verts = _sp.ndarray(network.num_throats(),dtype=object)
    for i in range(len(verts)):
        verts[i]=_sp.asarray(list(network["throat.vert_index"][i].values()))
    "Number of voxels in each direction"    
    lx = np.int((domain[0]/vox_len))
    ly = np.int((domain[1]/vox_len))
    lz = np.int((domain[2]/vox_len))
    "If the domain is less than 64000000 voxels we can generate the fibre image in one go"
    if lx*ly*lz <= 400**3:
        fibre_image = _get_fibre_image(network,geom_pores,vox_len,fibre_rad,verts)
    else:
        raise MemoryError("The network domain is too large to use the Voronoi Voxel method for pore volume")
    chunk_len = 100
    "If domain to big need to split into chunks to manage memory and multi-process"    
    if (lx > chunk_len) or (ly > chunk_len) or (lz > chunk_len):
        cx = np.ceil(lx/chunk_len).astype(int)
        cy = np.ceil(ly/chunk_len).astype(int)
        cz = np.ceil(lz/chunk_len).astype(int)
    else:
        cx = cy = cz = 1
    pore_chunks=[]
    for ci in range(cx):
        for cj in range(cy):
            for ck in range(cz):
                
                "Work out chunk range"
                cxmin = ci*chunk_len*vox_len
                cxmax = (ci+1)*chunk_len*vox_len
                cymin = cj*chunk_len*vox_len
                cymax = (cj+1)*chunk_len*vox_len
                czmin = ck*chunk_len*vox_len
                czmax = (ck+1)*chunk_len*vox_len
                crange=[cxmin,cxmax,cymin,cymax,czmin,czmax]
                "Find points within the range and get their throat vertices"
                cpores = network.pores()[(points[:,0]>=cxmin)*(points[:,0]<=cxmax)
                                        *(points[:,1]>=cymin)*(points[:,1]<=cymax)
                                        *(points[:,2]>=czmin)*(points[:,2]<=czmax)]
                pore_chunks.append([cpores,crange])
                #print("Num Pores in Chunk: "+str(len(cpores)))
                #print("Chunk Range: "+str(np.around(crange,4)))
                #print("Vertex Range: "+str(np.around(vo.vertex_dimension(network=network,face1=cpores,parm='minmax'),4)))
        
    #p = Pool(4)
    #chunk_vols = p.map(_get_voxel_volume, pore_chunks)
    for chunk_id in range(len(pore_chunks)):
        print("Processing Chunk: "+str(chunk_id+1)+" of "+str(len(pore_chunks)))
        chunk_pvols,chunk_fvols,chunk_pores = _get_voxel_volume(network,pore_chunks[chunk_id],vox_len,fibre_rad,verts,fibre_image) #Command to run when not in parallel
        volume[chunk_pores] +=chunk_pvols*voxel # this volume may not be the entire pore volume as some pores span multiple chunks, hence the addition
        fibre_vox[chunk_pores] +=chunk_fvols
        #volume[pore_chunks[chunk_id][1]] +=chunk_vols[chunk_id]
    geometry["pore.fibre_voxels"]=fibre_vox[geom_pores]
    return volume[geom_pores]

#if __name__ == '__main__':
#    freeze_support()