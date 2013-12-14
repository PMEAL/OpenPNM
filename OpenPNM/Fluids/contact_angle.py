
"""
module contact_angle
===============================================================================

"""
import scipy as sp

def constant(fluid,value=120,**params):
    fluid.pore_conditions['contact_angle'] = value

def na(fluid,**params):
    fluid.pore_conditions['contact_angle'] = -1





