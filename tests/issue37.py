# -*- coding: utf-8 -*-
"""
Created on Wed Oct 02 12:23:00 2013

@author: jeff
"""

import OpenPNM

pn = OpenPNM.Geometry.Cubic().generate()

#Now transate the network by 1 in all 3 dimensions
OpenPNM.Geometry.Cubic.translate_coordinates(pn,translation=[1,1,1])

