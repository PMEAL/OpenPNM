import sys
projpath =  __file__.rsplit('\\',2)[0]
sys.path.append(projpath)

import OpenPNM
import scipy as sp
from time import clock
start=clock()

#======================================================================
'''Generate Fluids'''
#======================================================================
#Define the fluids properties
air_recipe = {
'name': 'air',
'Pc': 3.771e6, #Pa
'Tc': 132.65,  #K
'MW': 0.0291,  #kg/mol
'diffusivity': {'method': 'Fuller',
                'MA': 0.03199,
                'MB': 0.0291,
                'vA': 16.3,
                'vB': 19.7},
'viscosity': {'method': 'Reynolds',
              'uo': 0.001,
              'b': 0.1},
'molar_density': {'method': 'ideal_gas',
                  'R': 8.314},
'surface_tension': {'method': 'constant',
                    'value': 0},
'contact_angle': {'method': 'na'},
}
water_recipe = {
'name': 'water',
'Pc': 2.206e6, #Pa
'Tc': 647,     #K
'MW': 0.0181,  #kg/mol
'diffusivity': {'method': 'constant',
                'value': 1e-12},
'viscosity': {'method': 'constant',
              'value': 0.001},
'molar_density': {'method': 'constant',
                  'value': 44445},
'surface_tension': {'method': 'Eotvos',
                    'k': 2.25e-4},
'contact_angle': {'method': 'constant',
                  'value': 120},
}
#Create fluids
air = OpenPNM.Fluids.GenericFluid(loglevel=50).create(air_recipe)
water= OpenPNM.Fluids.GenericFluid(loglevel=50).create(water_recipe)
#set water and air as a fluid pair
water.set_pair(air)
#Set desired base conditions in the Fluids
air.pore_conditions['temperature'] = 353
air.pore_conditions['pressure'] = 101325
water.pore_conditions['temperature'] = 353
water.pore_conditions['pressure'] = 101325
#Update Fluids to the new conditions
water.regenerate()
air.regenerate()


