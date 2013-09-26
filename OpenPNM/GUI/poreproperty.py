from __future__ import absolute_import
import OpenPNM
import numpy as np

def InletsAndOutlets(network, inlet_face=1, outlet_face=6):
    mask = network.pore_properties['type']==inlet_face
    inlets = network.pore_properties['numbering'][mask]
    network.pore_properties['inlets'][range(network.get_num_pores())]=0
    network.pore_properties['inlets'][inlets]=1
    mask = network.pore_properties['type']==outlet_face
    inlets = network.pore_properties['numbering'][mask]
    network.pore_properties['outlets'][range(network.get_num_pores())]=0
    network.pore_properties['outlets'][inlets]=1
    return {'network': network}
