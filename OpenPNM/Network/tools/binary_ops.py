# -*- coding: utf-8 -*-
"""
===============================================================================
Network.tools.binary_ops: Assorted topological manipulation methods
===============================================================================

"""
import scipy as _sp
from OpenPNM.Base import Controller, logging
ctrl = Controller()
logger = logging()

def add_pores(network,pore_coords=[],labels=[]):
    r'''
    Add individual pores to the simulation from a list of coordinates.

    Parameters
    ----------
    pore_coords : array_like
        The [X,Y,Z] coordinates of the pores to add 
    labels : string, or list of strings, optional
        A list of labels to apply to the new pores and throats

    '''

    logger.debug('Extending network')
    Np_old = network.num_pores()
    Np = Np_old + int(_sp.size(pore_coords)/3)
    #Adjust 'all' labels
    del network['pore.all']
    network['pore.all'] = _sp.ones((Np,),dtype=bool)
    #Add coords and conns
    if pore_coords != []:
        coords = _sp.vstack((network['pore.coords'],pore_coords))
        network['pore.coords'] = coords
    for item in network.keys():
        if item.split('.')[1] not in ['coords','all']:
            if item.split('.')[0] == 'pore':
                if network[item].dtype == bool:
                    temp = network[item]
                    network[item] = _sp.zeros((Np,),dtype=bool)
                    network[item][temp] = True
                elif network[item].dtype == object:
                    temp = network[item]
                    network[item] = _sp.ndarray((Np,),dtype=object)
                    network[item][_sp.arange(0,_sp.shape(temp)[0])] = temp
                else:
                    temp = network[item]
                    try:
                        network[item] = _sp.ones((Np,_sp.shape(temp)[1]),dtype=float)*_sp.nan
                    except:
                        network[item] = _sp.ones((Np,),dtype=float)*_sp.nan
                network[item][_sp.arange(0,_sp.shape(temp)[0])] = temp
    #Apply labels, if supplied
    if labels != []:
        #Convert labels to list if necessary
        if type(labels) is str:
            labels = [labels]
        for label in labels:
            #Remove pore or throat from label, if present
            label = label.split('.')[-1]
            if pore_coords != []:
                Ps = _sp.r_[Np_old:Np]
                if 'pore.'+label not in network.labels():
                    network['pore.'+label] = False
                network['pore.'+label][Ps] = True


    self._update_network()