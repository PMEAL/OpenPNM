# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 15:38:49 2014

@author: Jackie
"""

import scipy as sp
import matplotlib.pyplot as plt
import scipy.stats as spst

#generating my list of diameters
lmbda = 9e-6
k = 3.5
bmin = 9e-6
xmax = .9

Np = 100
x_values = [1.0*X/Np for X in range(1, Np)]
x_values2 = [1 - x_values[i] for i in range(len(x_values))]
values = [lmbda*(-sp.log(1-x_values[x]*xmax))**(-1/k) + bmin for x in range(len(x_values))]

#shape = 1
#loc = 6e-6
#scale = 20e-6
#lowestRMSE = 100
#bestFit = []
#
##preliminary narrowing
#for x in range(100):
#    
#    bestShape = shape
#    
#    #working on shape
#    for x in range(100, 500):
#        shape = x/100  
#        #generating other list of diameters
#        P = spst.weibull_min(shape,loc = loc,scale=scale)
#        values2 = P.ppf(x_values2)
#        RMSE = sp.sqrt(sp.mean((values-values2)**2))
#        if(RMSE < lowestRMSE):
#            lowestRMSE = RMSE
#            bestFit = values2
#            bestShape = shape
#        
#    shape = bestShape
#    
#    bestLoc = loc
#
#    #working on loc
#    for x in range(5, 250):
#        loc = x*1e-7 
#        #generating other list of diameters
#        P = spst.weibull_min(shape,loc = loc,scale = scale)
#        values2 = P.ppf(x_values2)
#        RMSE = sp.sqrt(sp.mean((values-values2)**2))
#        if(RMSE < lowestRMSE):
#            lowestRMSE = RMSE
#            bestFit = values2
#            bestLoc = loc
#        
#    loc = bestLoc
#
#    bestScale = scale
#    
#    #working on loc
#    for x in range(5, 250):
#        scale = x*1e-7  
#        #generating other list of diameters
#        P = spst.weibull_min(shape,loc = loc,scale = scale)
#        values2 = P.ppf(x_values2)
#        RMSE = sp.sqrt(sp.mean((values-values2)**2))
#        if(RMSE < lowestRMSE):
#            lowestRMSE = RMSE
#            bestFit = values2
#            bestScale = scale
#        
#    scale = bestScale
#
#print("lowest RMSE: ", lowestRMSE)  
#plt.plot(x_values, values, color = 'k', linestyle = '-')
#plt.plot(x_values, bestFit, color = 'r', linestyle = '-')
#plt.show()

#My results from preliminary narrowing
#temp_shape = 0.987
#temp_loc = 1.5799999999999998e-05
#temp_scale = 5.0999999999999995e-06

#P = spst.weibull_min(temp_shape,loc = temp_loc,scale=temp_scale)
#bestFit = P.ppf(x_values2)
#RMSE = sp.sqrt(sp.mean((values-bestFit)**2))
#
#My results from preliminary narrowing
#print("lowest RMSE: ", RMSE)  
#plt.plot(x_values, values, color = 'k', linestyle = '-')
#plt.plot(x_values, bestFit, color = 'r', linestyle = '-')
#plt.show()

#secondary narrowing
#lowestRMSE = 6.41511262453e-07
#scale = temp_scale
#loc = temp_loc
#shape = temp_shape
#
#for x in range(5):
#    bestShape = shape
#    
#    #working on shape
#    for x in range(-500, 500):
#        shape = temp_shape + x/1000  
#        #generating other list of diameters
#        P = spst.weibull_min(shape,loc = loc,scale=scale)
#        values2 = P.ppf(x_values2)
#        RMSE = sp.sqrt(sp.mean((values-values2)**2))
#        if(RMSE < lowestRMSE):
#            lowestRMSE = RMSE
#            bestFit = values2
#            bestShape = shape
#        
#    shape = bestShape
#    
#    bestLoc = loc
#
#    #working on loc
#    for x in range(-500, 500):
#        loc = temp_loc + x/1000 
#        #generating other list of diameters
#        P = spst.weibull_min(shape,loc = loc,scale = scale)
#        values2 = P.ppf(x_values2)
#        RMSE = sp.sqrt(sp.mean((values-values2)**2))
#        if(RMSE < lowestRMSE):
#            lowestRMSE = RMSE
#            bestFit = values2
#            bestLoc = loc
#        
#    loc = bestLoc
#
#    bestScale = scale
#    
#    #working on loc
#    for x in range(-500, 500):
#        scale = temp_scale + x/1000 
#        #generating other list of diameters
#        P = spst.weibull_min(shape,loc = loc,scale = scale)
#        values2 = P.ppf(x_values2)
#        RMSE = sp.sqrt(sp.mean((values-values2)**2))
#        if(RMSE < lowestRMSE):
#            lowestRMSE = RMSE
#            bestFit = values2
#            bestScale = scale
#        
#    scale = bestScale
#
#print("lowest RMSE: ", lowestRMSE)  
#plt.plot(x_values, values, color = 'k', linestyle = '-')
#plt.plot(x_values, bestFit, color = 'r', linestyle = '-')
#plt.show()

#final results
P = spst.weibull_min(.987,loc = 1.58e-5,scale=5.01e-6)
values2 = P.ppf(x_values2)
RMSE = sp.sqrt(sp.mean((values-values2)**2))
print("RMSE: ", RMSE)
p1, = plt.plot(x_values, values, color = 'r', linestyle = '-')
p2, = plt.plot(x_values, values2, color = 'b', linestyle = '-')
plt.legend([p1, p2],
           [r'$\lambda[-\ln(1-xx_{max})]^{-1/k} + b_{min} $' + '\n' + r'$\lambda = 9*10^{-9} $' + '\n' + r'$x_{max} = .9$' + '\n' + r'$k = 3.5$' + '\n' + r'$b_{min} = 9*10^{-6}$', 
           r'$scipy\hspace{1}weibull\hspace{1}min$' +'\n' + r'$RMSE: 6.31616833715*10^{-7}$' + '\n' + r'$shape = .987$' + '\n' + r'$loc = 1.58*10^{-5}$'+ '\n' + r'$scale = 5.01*10^{-6}$'], 
            loc = 'upper right') 
plt.figtext(.48, .3, 'Note: x values for weibull are really 1-x', bbox=dict(facecolor='white', alpha=0.5))
plt.xlabel("random x between 0 and 1")
plt.ylabel("pore diameter")
plt.show()




