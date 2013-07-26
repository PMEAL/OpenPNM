# -*- coding: utf-8 -*-
"""
Created on Tue Jul 09 08:32:53 2013

@author: Jeff
"""

import scipy.spatial as sptl
import scipy.sparse as sprs
Lx = params['domain_size'][0]
Ly = params['domain_size'][1]
Lz = params['domain_size'][2]
Np = pn.get_num_pores()
btype = params['btype']

#Translate internal pores to each face of domain
pcoords0 = pn.pore_properties['coords']
pcoords1 = pcoords0 - [0,0,Lz]
pcoords2 = pcoords0 - [Lx,0,0]    
pcoords3 = pcoords0 - [0,Ly,0]
pcoords4 = pcoords0 + [0,Ly,0]
pcoords5 = pcoords0 + [Lx,0,0]
pcoords6 = pcoords0 + [0,0,Lz]
pts = np.vstack((pcoords0,pcoords1,pcoords2,pcoords3,pcoords4,pcoords5,pcoords6))

#Use some twisted logic to get bval list of 1 for boundary and -1 for periodic faces
bval = [0, 1, 2, 3, 4, 5, 6]*(np.array([0, btype[2], btype[0], btype[1], btype[1], btype[0], btype[2]])*-2+1)
ptype0 = np.zeros((Np,),dtype=int)
ptype1 = np.ones_like(ptype0)*bval[1]
ptype2 = np.ones_like(ptype0)*bval[2]
ptype3 = np.ones_like(ptype0)*bval[3]
ptype4 = np.ones_like(ptype0)*bval[4]
ptype5 = np.ones_like(ptype0)*bval[5]
ptype6 = np.ones_like(ptype0)*bval[6]
#typ is a list of pore types (+2 & +5 = regular x face, -1 & -6 = periodic z face)
typ = np.hstack((ptype0,ptype1,ptype2,ptype3,ptype4,ptype5,ptype6))

pnum0 = pn.pore_properties['numbering'][pn.pore_properties['type']==0]
#num contains the internal ID number of the boundary pores (for connecting periodic points)
num = np.tile(pnum0,7)

Tri = sptl.Delaunay(pts)
adjmat = sprs.lil_matrix((np.shape(pts)[0],np.shape(pts)[0]),dtype=int)
for i in np.arange(0,np.shape(Tri.simplices)[0]):
    #Keep only simplices that are fully in real domain
    adjmat[Tri.simplices[i],Tri.simplices[i]] = 1
adjmat = sprs.triu(adjmat,k=1,format="lil")
for i in np.arange(0,Np):
    #Add periodic throats to the netowrk (if any)
    tpore2 = num[adjmat.rows[i]][typ[adjmat.rows[i]]<0]
    tpore1 = np.ones_like(tpore2,dtype=int)*i
    pn.throat_properties['connections'] = np.concatenate((pn.throat_properties['connections'],np.vstack((tpore1,tpore2)).T),axis=0)
    pn.throat_properties['type'] = np.concatenate((pn.throat_properties['type'],typ[adjmat.rows[i]][typ[adjmat.rows[i]]<0]))
    #Add boundary pores and throats to the network
    newporetyps = np.unique(typ[adjmat.rows[i]][typ[adjmat.rows[i]]>0])
    newporenums = np.r_[pn.get_num_pores():pn.get_num_pores()+np.size(newporetyps)]
    tpore2 = newporenums
    tpore1 = np.ones_like(tpore2,dtype=int)*i
    pn.throat_properties['connections'] = np.concatenate((pn.throat_properties['connections'],np.vstack((tpore1,tpore2)).T),axis=0)
    pn.throat_properties['type'] = np.concatenate((pn.throat_properties['type'],newporetyps),axis=0)            
    pn.pore_properties['type'] = np.concatenate((pn.pore_properties['type'],newporetyps),axis=0)
    bcoords = np.zeros((7,3),dtype=float)
    bcoords[1,:] = [pn.pore_properties['coords'][i,0], pn.pore_properties['coords'][i,1], 0]
    bcoords[2,:] = [0, pn.pore_properties['coords'][i,1], pn.pore_properties['coords'][i,2]]
    bcoords[3,:] = [pn.pore_properties['coords'][i,0], 0, pn.pore_properties['coords'][i,2]]
    bcoords[4,:] = [pn.pore_properties['coords'][i,0], pn.pore_properties['coords'][i,1], Lz]
    bcoords[5,:] = [Lx, pn.pore_properties['coords'][i,1], pn.pore_properties['coords'][i,2]]
    bcoords[6,:] = [pn.pore_properties['coords'][i,0], Ly, pn.pore_properties['coords'][i,2]]
    newporecoords = bcoords[newporetyps,:]
    pn.pore_properties['coords'] = np.concatenate((pn.pore_properties['coords'],newporecoords),axis=0)
#Reset number of pores and throats (easier than tracking it)
pn.pore_properties['numbering'] = np.r_[0:np.size(pn.pore_properties['type'])]
pn.throat_properties['numbering'] = np.r_[0:np.size(pn.throat_properties['type'])]