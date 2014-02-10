import sys, os
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if sys.path[1] != parent_dir:
    sys.path.insert(1, parent_dir)

import doctest
verbose = False
report = True

#------------------------------------------------------------------------------
'''Test Base Classes'''
#------------------------------------------------------------------------------
import OpenPNM.Base.__OpenPNMbase__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

import OpenPNM.Base.__Tools__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

#------------------------------------------------------------------------------
'''Test Network Classes'''
#------------------------------------------------------------------------------
import OpenPNM.Network.__GenericNetwork__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

import OpenPNM.Network.__Cubic__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

import OpenPNM.Network.__Delaunay__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

import OpenPNM.Network.__Template__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

#------------------------------------------------------------------------------
'''Test Geometry Classes'''
#------------------------------------------------------------------------------
import OpenPNM.Geometry.__GenericGeometry__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

import OpenPNM.Geometry.__StickBall__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

#------------------------------------------------------------------------------
'''Test Fluid Classes'''
#------------------------------------------------------------------------------
import OpenPNM.Fluids.__GenericFluid__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

import OpenPNM.Fluids.__Air__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)

import OpenPNM.Fluids.__Water__ as mod_to_test
doctest.testmod(mod_to_test,report=report,verbose=verbose)
