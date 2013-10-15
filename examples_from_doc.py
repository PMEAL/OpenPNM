# -*- coding: utf-8 -*-

# Getting started

import OpenPNM
gn1 = OpenPNM.Geometry.Cubic()
pn1 = gn1.generate()
print pn1
pn2 = OpenPNM.Geometry.Cubic().generate()

#Network Architecture
pn = OpenPNM.Geometry.Cubic().generate(divisions=[3,3,3],lattice_spacing=[1])
print pn
pn.print_overview()
print pn.pore_properties.keys()
pn.pore_conditions['temperature'] = 80.0
pn.print_overview()