# -*- coding: utf-8 -*-

# Getting started

import OpenPNM
gn1 = OpenPNM.Geometry.Cubic()
pn1 = gn1.generate()
print pn1.pore_properties['diameter'][0]
print pn1
pn2 = OpenPNM.Geometry.Cubic().generate()
print pn1.pore_properties['diameter'][0]
print pn2.pore_properties['diameter'][0]


#Network Architecture
pn = OpenPNM.Geometry.Cubic().generate(divisions=[3,3,3],lattice_spacing=[1])
pn.print_overview()