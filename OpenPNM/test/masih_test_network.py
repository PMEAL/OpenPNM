# -*- coding: utf-8 -*-
"""
Created on Wed Mar 06 16:00:41 2013

@author: Masih
"""

import numpy as np
import matplotlib.pyplot as plt
kappa_p = 6
lamda_p = 9

Nx = 30
Ny = 30
Np = Nx*Ny
Nt=(Nx-1)*(Ny-1)*2+(Nx-1)+(Ny-1)
#Assign pore size to each pore
psize = np.random.rand(Np)
#Apply Weibull distribution
pdia = 2*lamda_p*(-1*np.log(1.0-(psize*0.9)))**(1.0/kappa_p)
#Limit maximum pore size to 25        
pdia[pdia>25] = 25
#Limit minimum pore size to 15 
pdia[pdia<15] = 15
kappa_t = 4
lamda_t = 3
tsize = np.random.rand(Nt)
tdia = 2*lamda_t*(-1*np.log(1.0-(tsize*0.9)))**(1.0/kappa_t)
#Limit maximum throat size to 25        
tdia[tdia>10] = 10
#Limit minimum throat size to 15 
tdia[tdia<5] = 5       
#Initialize tpore1 and tpore larger than needed, the trim later
tpore1 = np.zeros(Nt,int)
tpore2 = np.zeros(Nt,int)
xcoord=np.zeros(Np,int)
ycoord=np.zeros(Np,int)
counter=0
counti=0
for j in range(Ny-1):
    for k in range(Nx-1):
        xcoord[counter]=k
        ycoord[counter]=j
        porenum=Ny*j+k
        tpore1[counti] =porenum
        tpore2[counti]=porenum+1
        counti=counti+1
        tpore1[counti]=porenum
        tpore2[counti]=porenum+Ny
        counti=counti+1
        counter=counter+1
    k=Nx-1
    xcoord[counter]=k
    ycoord[counter]=j
    porenum=Ny*j+k
    tpore1[counti]=porenum
    tpore2[counti]=porenum+Ny
    counti=counti+1
    counter=counter+1
j=Ny-1
for k in range(Nx-1):
    xcoord[counter]=k
    ycoord[counter]=j
    porenum=Ny*j+k
    tpore1[counti] =porenum
    tpore2[counti]=porenum+1
    counti=counti+1
    counter=counter+1
    
xcoord[899]=29
ycoord[899]=29
#         tpore1
#throut_r=np.random.randint(5,10, [Nx, Ny])
#pore_r=np.random.randint(15,25,[Nx, Ny])
##plt.imshow(pore_r,interpolation='nearest')    
##plt.show()



import OpenPNM.NET as NET
import scipy as sp
pn = NET.BaseNet(num_pores = Np, num_throats = sp.size(tpore1))
pn.set_throat_property(name = 'connections',ndarray = sp.transpose(sp.array([tpore1,tpore2])))
pn.set_throat_property(name= 'dia',ndarray=tdia)
pn.set_throat_property(name = 'pore1', ndarray=tpore1)
pn.set_throat_property(name = 'pore2', ndarray=tpore2)
pn.set_pore_property(name = 'dia', ndarray = pdia)
pn.set_pore_property(name = 'coordx', ndarray=xcoord)
pn.set_pore_property(name = 'coordy', ndarray=ycoord)
pn.set_pore_property(name = 'coordz', ndarray=ycoord*0)

import OpenPNM.Visualization as Vis
Vis.NetToVtp(net = pn)

