import numpy as np

def pnm_2_salome(sphere_c, sphere_r, cylinder_head, cylinder_tail, cylinder_r):
    n = len(sphere_r)
    m = len(cylinder_r)
    vertex_id = 0
    
    fuse_list = np.array([])
    
    for i in range(n):
        vertex_s = geompy.MakeVertex(sphere_c[i])
        sphere = geompy.MakeSpherePntR(vertex_s, sphere_r[i])
        
        geompy.addToStudy(vertex_s, "Vertex_{}".format(vertex_id))
        vertex_id += 1
        geompy.addToStudy(sphere, "Sphere_{}".format(i))
        
        fuse_list = np.append(fuse_list, sphere)

    for i in range(m):
        vertex_c1 = geompy.MakeVertex(cylinder_head[i])
        vertex_c2 = geompy.MakeVertex(cylinder_tail[i])
        line = geompy.MakeLineTwoPnt(vertex_c1, vertex_c2)
        length = abs(np.linalg.norm(cylinder_head[i]) -
                     np.linalg.norm(cylinder_tail[i]))
        cylinder = geompy.MakeCylinder(vertex_c1, line, cylinder_r[i], length)
        
        geompy.addToStudy(vertex_c1, "Vertex_{}".format(vertex_id))
        vertex_id += 1
        geompy.addToStudy(vertex_c2, "Vertex_{}".format(vertex_id))
        vertex_id += 1
        geompy.addToStudy(line, "Line_{}".format(i))
        geompy.addToStudy(cylinder, "Cylinder_{}".format(i))
        
        fuse_list = np.append(fuse_list, cylinder)
    
    
    geompy.addToStudy(fuse_list, "Fuse")

