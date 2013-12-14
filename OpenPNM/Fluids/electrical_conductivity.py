
"""
module electrical_conductance
===============================================================================

"""
import scipy as sp

def constant(fluid,value=1,**params):
    fluid.pore_conditions['electrical_conductance'] = value

def na(fluid,**params):
    fluid.pore_conditions['electrical_conductance'] = -1

