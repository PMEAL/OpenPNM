# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 09:07:48 2021

@author: xu kai
"""

#=================================
#prepare

import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
import scipy as sp
import openpnm.models as mods
import openpnm.io.VTK as iovtk
from bimodal_distribution import bimodal_distribution
from find_organic_throats import find_organic_throats

# network model
pn = op.network.Cubic(shape=[10, 10, 10], spacing=0.0001)

# geometry model
geo = op.geometry.GenericGeometry(network=pn, pores=pn.pores(),
                                  throats=pn.throats(), name='geo')

geo.add_model(propname='pore.seed',
             model=mods.misc.random,
             element='pore',
             num_range=[0.2, 0.7],
             seed=None)

geo.add_model(propname='pore.max_size',
              model=mods.geometry.pore_size.largest_sphere,
              iters=10)

geo.add_model(propname='pore.diameter',
              model=mods.misc.product,
              prop1='pore.max_size',
              prop2='pore.seed')

# generate a bimodal distribution
geo['pore.diameter'] = bimodal_distribution()

geo.add_model(propname='throat.diameter',
              model=mods.geometry.throat_size.from_neighbor_pores,
              pore_prop='pore.diameter',
              mode='min')

geo.add_model(propname='pore.area',
              model=mods.geometry.pore_area.sphere)

geo.add_model(propname='pore.volume',
              model=mods.geometry.pore_volume.sphere)

geo.add_model(propname='throat.endpoints',
                       model=mods.geometry.throat_endpoints.spherical_pores,
                       pore_diameter='pore.diameter',
                       throat_diameter='throat.diameter')

geo.add_model(propname='throat.conduit_lengths',
                       model=mods.geometry.throat_length.conduit_lengths,
                       throat_endpoints='throat.endpoints',
                       throat_length='throat.length')

geo.add_model(propname='throat.length',
              model=mods.geometry.throat_length.piecewise)

geo.add_model(propname='throat.volume',
              model=mods.geometry.throat_volume.cylinder)

geo.add_model(propname='throat.area',
              model=mods.geometry.throat_area.cylinder)

geo.add_model(propname='throat.surface_area',
              model=mods.geometry.throat_surface_area.cylinder)

geo['throat.radius'] = geo['throat.diameter'] / 2

# shape factor
geo.add_model(propname='throat.perimeter',
              model=mods.geometry.throat_perimeter.cylinder,
              throat_diameter='throat.diameter')

geo.add_model(propname='throat.area',
              model=mods.geometry.throat_area.cylinder)

geo.add_model(propname='throat.shape_factor',
              model=mods.geometry.throat_capillary_shape_factor.mason_morrow,
              throat_perimeter='throat.perimeter',
              throat_area='throat.area')

# plan a for shape factor
geo['pore.perimeter'] = np.pi * geo['pore.diameter']
geo['pore.shape_factor'] = geo['pore.area'] / geo['pore.perimeter'] ** 2

# plan b for shape factor
#geo['pore.perimeter'] = geo['pore.volume'] / geo['pore.area']**(3/2)

# organic pores & inorganic pores
# only need to identify organic pores, than the left pores are inorganic pores
geo['pore.organic'] = geo['pore.diameter'] <= 50e-9

organic_pores = [pore_index for pore_index in geo.pores()
                if geo['pore.organic'][pore_index]==True]
organic_pores = np.array(organic_pores)

# find organic throats
organic_throats = find_organic_throats(organic_pores,
                                       geo['throat.conns'], geo.Nt)

# physics
water = op.phases.Water(network=pn)

# viscosity settings
alpha_o = 1.1  # viscosity ratio organic pores = \mu_a / \mu, range(1-2.5)
alpha_i = 0.9  # viscosity ratio in inorganic pores, range(0.5-1)

# initialize
water['pore.viscosity'] = 3.6e-3
water['throat.viscosity'] = 3.6e-3
water['pore.viscosity_a'] = water['pore.viscosity'] * alpha_i
water['throat.viscosity_a'] = water['throat.viscosity'] * alpha_i

water['pore.viscosity_a'][organic_pores] =\
                water['pore.viscosity'][organic_pores] * alpha_o
water['throat.viscosity_a'][organic_throats[1]] =\
                water['throat.viscosity'][organic_throats[1]] * alpha_o

# slip length
Ls_o = 60e-9 # organic slip length, range 0-250 nm
Ls_i = 50e-9 # organic slip length, range 0-60 nm

# slip length of inorganic pores and throats
water['pore.l_sd'] = Ls_i /  np.sqrt(geo['pore.area'])
water['throat.l_sd'] = Ls_i /  np.sqrt(geo['throat.area'])

# # dimensionless slip length of organic pores and throats
water['pore.l_sd'][organic_pores] =\
                Ls_o / np.sqrt(geo['pore.area'][organic_pores])
water['throat.l_sd'][organic_throats[1]] =\
                Ls_o /  np.sqrt(geo['throat.area'][organic_throats[1]])

# thickness of adsorption layer
ha = 1.8e-9
geo['pore.area_a'] = geo['pore.area'] - \
                     np.pi * (geo['pore.diameter'] / 2 - ha) ** 2
geo['throat.area_a'] = geo['throat.area'] - \
                       np.pi * (geo['throat.radius'] - ha) ** 2

# hydraulic conductance


# porosity
porosity = (geo['pore.volume'].sum() + geo['throat.volume'].sum()) / 0.001**3















