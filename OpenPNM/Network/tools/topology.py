# -*- coding: utf-8 -*-
"""
===============================================================================
Network.tools.topology: Assorted topological manipulation methods
===============================================================================

"""
import scipy as _sp
from OpenPNM.Base import Controller, logging
ctrl = Controller()
logger = logging.getLogger(__name__)

def extend(network,pore_coords=[],throat_conns=[],labels=[]):
    r'''
    Add individual pores (or throats) to the network from a list of coords
    or conns.

    Parameters
    ----------
    pore_coords : array_like
        The coordinates of the pores to add
    throat_conns : array_like
        The throat connections to add
    labels : string, or list of strings, optional
        A list of labels to apply to the new pores and throats

    Notes
    -----
    This needs to be enhanced so that it increases the size of all pore
    and throat props and labels on ALL associated objects.  At the moment
    if throws an error is there are ANY associated objects.

    '''
    if (network._geometries != []):
        raise Exception('Network has active Geometries, cannot proceed')
    if (network._phases != []):
        raise Exception('Network has active Phases, cannot proceed')

    logger.debug('Extending network')
    Np_old = network.num_pores()
    Nt_old = network.num_throats()
    Np = Np_old + int(_sp.size(pore_coords)/3)
    Nt = Nt_old + int(_sp.size(throat_conns)/2)
    #Adjust 'all' labels
    del network['pore.all'], network['throat.all']
    network['pore.all'] = _sp.ones((Np,),dtype=bool)
    network['throat.all'] = _sp.ones((Nt,),dtype=bool)
    #Add coords and conns
    if pore_coords != []:
        coords = _sp.vstack((network['pore.coords'],pore_coords))
        network['pore.coords'] = coords
    if throat_conns != []:
        conns = _sp.vstack((network['throat.conns'],throat_conns))
        network['throat.conns'] = conns
    for item in network.keys():
        if item.split('.')[1] not in ['coords','conns','all']:
            if item.split('.')[0] == 'pore':
                N = Np
            else:
                N = Nt
            if network[item].dtype == bool:
                temp = network[item]
                network[item] = _sp.zeros((N,),dtype=bool)
                network[item][temp] = True
            elif network[item].dtype == object:
                temp = network[item]
                network[item] = _sp.ndarray((N,),dtype=object)
                network[item][_sp.arange(0,_sp.shape(temp)[0])] = temp
            else:
                temp = network[item]
                try:
                    network[item] = _sp.ones((N,_sp.shape(temp)[1]),dtype=float)*_sp.nan
                except:
                    network[item] = _sp.ones((N,),dtype=float)*_sp.nan
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
            if throat_conns != []:
                Ts = _sp.r_[Nt_old:Nt]
                if 'throat.'+label not in network.labels():
                    network['throat.'+label] = False
                network['throat.'+label][Ts] = True

    network._update_network()
    
def trim(self, pores=[], throats=[]):
    '''
    Remove pores (or throats) from the network.

    Parameters
    ----------
    pores (or throats) : array_like
        A boolean mask of length Np (or Nt) or a list of indices of the
        pores (or throats) to be removed.

    Notes
    -----
    Trimming only adjusts Phase, Geometry, and Physics objects. Trimming a
    Network that has already been used to run simulations will break those
    simulation objects.

    Examples
    --------
    >>> import OpenPNM
    >>> pn = OpenPNM.Network.TestNet()
    >>> pn.Np
    125
    >>> pn.Nt
    300
    >>> pn.trim(pores=[1])
    >>> pn.Np
    124
    >>> pn.Nt
    296

    '''
    for net in self.controller.networks():
        if net._parent is self:
            raise Exception('This Network has been cloned, cannot trim')

    if pores != []:
        pores = _sp.array(pores,ndmin=1)
        Pkeep = _sp.ones((self.num_pores(),),dtype=bool)
        Pkeep[pores] = False
        Tkeep = _sp.ones((self.num_throats(),),dtype=bool)
        Ts = self.find_neighbor_throats(pores)
        if len(Ts)>0:
            Tkeep[Ts] = False
    elif throats != []:
        throats = _sp.array(throats,ndmin=1)
        Tkeep = _sp.ones((self.num_throats(),),dtype=bool)
        Tkeep[throats] = False
        Pkeep = self['pore.all'].copy()
    else:
        logger.warning('No pores or throats recieved')
        return

    # Trim all associated objects
    for item in self._geometries+self._physics+self._phases:
        Pnet = self['pore.'+item.name]*Pkeep
        Tnet = self['throat.'+item.name]*Tkeep
        temp = self.map_pores(pores=_sp.where(Pnet)[0],target=item,return_mapping=True)
        Ps = temp['target']
        temp = self.map_throats(throats=_sp.where(Tnet)[0],target=item,return_mapping=True)
        Ts = temp['target']
        # Then resize 'all
        item.update({'pore.all' : _sp.ones((_sp.sum(Pnet),),dtype=bool)})
        item.update({'throat.all' : _sp.ones((_sp.sum(Tnet),),dtype=bool)})
        # Overwrite remaining data and info
        for key in list(item.keys()):
            if key.split('.')[1] not in ['all']:
                temp = item.pop(key)
                if key.split('.')[0] == 'throat':
                    logger.debug('Trimming {a} from {b}'.format(a=key,b=item.name))
                    item[key] = temp[Ts]
                if key.split('.')[0] == 'pore':
                    logger.debug('Trimming {a} from {b}'.format(a=key,b=item.name))
                    item[key] = temp[Ps]

    #Remap throat connections
    Pmap = _sp.ones((self.Np,),dtype=int)*-1
    Pmap[Pkeep] = _sp.arange(0,_sp.sum(Pkeep))
    tpore1 = self['throat.conns'][:,0]
    tpore2 = self['throat.conns'][:,1]
    Tnew1 = Pmap[tpore1[Tkeep]]
    Tnew2 = Pmap[tpore2[Tkeep]]
    #Write 'all' label specifically
    self.update({'throat.all' : _sp.ones((_sp.sum(Tkeep),),dtype=bool)})
    self.update({'pore.all' : _sp.ones((_sp.sum(Pkeep),),dtype=bool)})
    # Write throat connections specifically
    self.update({'throat.conns' : _sp.vstack((Tnew1,Tnew2)).T})
    # Overwrite remaining data and info
    for item in list(self.keys()):
        if item.split('.')[-1] not in ['conns','all']:
            temp = self.pop(item)
            if item.split('.')[0] == 'throat':
                logger.debug('Trimming {a} from {b}'.format(a=item,b=self.name))
                self[item] = temp[Tkeep]
            if item.split('.')[0] == 'pore':
                logger.debug('Trimming {a} from {b}'.format(a=item,b=self.name))
                self[item] = temp[Pkeep]

    #Reset network graphs
    self._update_network(mode='regenerate')

    #Check Network health
    health = self.check_network_health()
    if health['trim_pores'] != []:
        logger.warning('Isolated pores exist!  Run check_network_health to ID which pores to remove.')
        pass