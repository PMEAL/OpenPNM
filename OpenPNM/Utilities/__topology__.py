# -*- coding: utf-8 -*-
"""
===============================================================================
Network.tools.topology: Assorted topological manipulation methods
===============================================================================

"""
import scipy as _sp
import scipy.sparse as _sprs
import scipy.spatial as _sptl
from OpenPNM.Base import logging as _logging
from OpenPNM.Base import Controller as _controller
logger = _logging.getLogger(__name__)
_ctrl = _controller()

class topology(object):

    def extend(self,network,pore_coords=[],throat_conns=[],labels=[]):
        r'''
        Add individual pores and/or throats to the network from a list of coords
        or conns.  This is an in-place operation, meaning the received Network
        object will be altered directly.
    
        Parameters
        ----------
        network : OpenPNM Network Object
            The Network to which pores or throats should be added
        pore_coords : array_like
            The coordinates of the pores to add
        throat_conns : array_like
            The throat connections to add
        labels : string, or list of strings, optional
            A list of labels to apply to the new pores and throats
    
        Notes
        -----
        This needs to be enhanced so that it increases the size of all pore
        and throat props and labels on ALL associated Phase objects.  At the
        moment it throws an error is there are any associated Phases.
    
        '''
        if (network._phases != []):
            raise Exception('Network has active Phases, cannot proceed')
    
        logger.info('Extending network')
        Np_old = network.num_pores()
        Nt_old = network.num_throats()
        Np = Np_old + int(_sp.size(pore_coords)/3)
        Nt = Nt_old + int(_sp.size(throat_conns)/2)
        # Adjust 'all' labels
        del network['pore.all'], network['throat.all']
        network['pore.all'] = _sp.ones((Np,),dtype=bool)
        network['throat.all'] = _sp.ones((Nt,),dtype=bool)
        # Add coords and conns
        if pore_coords != []:
            coords = _sp.vstack((network['pore.coords'],pore_coords))
            network['pore.coords'] = coords
        if throat_conns != []:
            conns = _sp.vstack((network['throat.conns'],throat_conns))
            network['throat.conns'] = conns
        # Increase size of any prop or label arrays on Network
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
        # Regnerate the adjacency matrices
        network._update_network()
        
    def trim(self,network, pores=[], throats=[]):
        '''
        Remove pores or throats from the network.  This is an in-place operation, 
        meaning the received Network object will be altered directly.
    
        Parameters
        ----------
        network : OpenPNM Network Object
            The Network from which pores or throats should be removed
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
        ctrl = network.controller
        for net in ctrl.networks():
            if net._parent is network:
                raise Exception('This Network has been cloned, cannot trim')
        if (pores != []) and (throats != []):
            raise Exception('Cannot delete pores and throats simultaneously')
        elif pores != []:
            pores = _sp.array(pores,ndmin=1)
            Pkeep = _sp.ones((network.num_pores(),),dtype=bool)
            Pkeep[pores] = False
            Tkeep = _sp.ones((network.num_throats(),),dtype=bool)
            Ts = network.find_neighbor_throats(pores)
            if len(Ts)>0:
                Tkeep[Ts] = False
        elif throats != []:
            throats = _sp.array(throats,ndmin=1)
            Tkeep = _sp.ones((network.num_throats(),),dtype=bool)
            Tkeep[throats] = False
            Pkeep = network['pore.all'].copy()
        else:
            logger.warning('No pores or throats recieved')
            return
    
        # Trim all associated objects
        for item in network._geometries+network._physics+network._phases:
            Pnet = network['pore.'+item.name]*Pkeep
            Tnet = network['throat.'+item.name]*Tkeep
            temp = network.map_pores(pores=_sp.where(Pnet)[0],target=item,return_mapping=True)
            Ps = temp['target']
            temp = network.map_throats(throats=_sp.where(Tnet)[0],target=item,return_mapping=True)
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
        Pmap = _sp.ones((network.Np,),dtype=int)*-1
        Pmap[Pkeep] = _sp.arange(0,_sp.sum(Pkeep))
        tpore1 = network['throat.conns'][:,0]
        tpore2 = network['throat.conns'][:,1]
        Tnew1 = Pmap[tpore1[Tkeep]]
        Tnew2 = Pmap[tpore2[Tkeep]]
        #Write 'all' label specifically
        network.update({'throat.all' : _sp.ones((_sp.sum(Tkeep),),dtype=bool)})
        network.update({'pore.all' : _sp.ones((_sp.sum(Pkeep),),dtype=bool)})
        # Write throat connections specifically
        network.update({'throat.conns' : _sp.vstack((Tnew1,Tnew2)).T})
        # Overwrite remaining data and info
        for item in list(network.keys()):
            if item.split('.')[-1] not in ['conns','all']:
                temp = network.pop(item)
                if item.split('.')[0] == 'throat':
                    logger.debug('Trimming {a} from {b}'.format(a=item,b=network.name))
                    network[item] = temp[Tkeep]
                if item.split('.')[0] == 'pore':
                    logger.debug('Trimming {a} from {b}'.format(a=item,b=network.name))
                    network[item] = temp[Pkeep]
    
        #Reset network graphs
        network._update_network(mode='regenerate')
    
        #Check Network health
        health = network.check_network_health()
        if health['trim_pores'] != []:
            logger.warning('Isolated pores exist!  Run check_network_health to ID which pores to remove.')
            pass
        
    def clone_pores(self,network,pores,apply_label=['clone'],mode='parents'):
        r'''
        Clones the specified pores and adds them to the network
    
        Parameters
        ----------
        network : OpenPNM Network Object
            The Network object to which the new pores are to be added
        pores : array_like
            List of pores to clone
        apply_labels : string, or list of strings
            The labels to apply to the clones, default is 'clone'
        mode : string
            Controls the connections between parents and clones.  Options are:
    
            - 'parents': (Default) Each clone is connected only to its parent
            - 'siblings': Clones are only connected to each other in the same manner as parents were connected
            - 'isolated': No connections between parents or siblings
        '''
        if (network._geometries != []):
            logger.warning('Network has active Geometries, new pores must be assigned a Geometry')
        if (network._phases != []):
            raise Exception('Network has active Phases, cannot proceed')
    
        logger.debug('Cloning pores')
        apply_label = list(apply_label)
        #Clone pores
        Np = network.num_pores()
        Nt = network.num_throats()
        parents = _sp.array(pores,ndmin=1)
        pcurrent = network['pore.coords']
        pclone = pcurrent[pores,:]
        pnew = _sp.concatenate((pcurrent,pclone),axis=0)
        Npnew = _sp.shape(pnew)[0]
        clones = _sp.arange(Np,Npnew)
        #Add clone labels to network
        for item in apply_label:
            if ('pore.'+item) not in network.keys():
                network['pore.'+item] = False
            if ('throat.'+item) not in network.keys():
                network['throat.'+item] = False
        #Add connections between parents and clones
        if mode == 'parents':
            tclone = _sp.vstack((parents,clones)).T
            self.extend(network=network,pore_coords=pclone,throat_conns=tclone)
        if mode == 'siblings':
            ts = network.find_neighbor_throats(pores=pores,mode='intersection')
            tclone = network['throat.conns'][ts] + network.num_pores()
            self.extend(network=network,pore_coords=pclone,throat_conns=tclone)
        if mode == 'isolated':
            self.extend(network=network,pore_coords=pclone)
        #Apply provided labels to cloned pores
        for item in apply_label:
            network['pore.'+item][network.pores('all')>=Np] = True
            network['throat.'+item][network.throats('all')>=Nt] = True
    
        # Any existing adjacency and incidence matrices will be invalid
        network._update_network()
        
    def stitch(self,network,donor,P_network,P_donor,method='delaunay',len_max=_sp.inf,label_suffix=''):
        r'''
        Stitches a second a network to the current network.
    
        Parameters
        ----------
        networK : OpenPNM Network Object
            The Network that will to which to donor Network will be attached
            
        donor : OpenPNM Network Object
            The Network to stitch on to the current Network
    
        pores_1 : array_like
            The pores on the current Network
    
        pores_2 : array_like
            The pores on the donor Network
            
        label_suffix : string or None
            Some text to append to each label in the donor Network before
            inserting them into the recipient.  The default is to append no 
            text, but a common option would be to append the donor Network's 
            name. To insert none of the donor labels, use None.
    
        len_max : float
            Set a length limit on length of new throats
    
        method : string (default = 'delaunay')
            The method to use when making pore to pore connections. Options are:
    
            - 'delaunay' : Use a Delaunay tessellation
            - 'nearest' : Connects each pore on the receptor network to its nearest pore on the donor network
            
        Notes
        -----
        Before stitching it is necessary to translate the pore coordinates of 
        one of the Networks so that it is positioned correctly relative to the
        other.  
        
        Examples
        --------
        >>> import OpenPNM
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = OpenPNM.Network.TestNet()
        >>> [pn.Np, pn.Nt]
        [125, 300]
        >>> [pn2.Np, pn2.Nt]
        [125, 300]
        >>> pn2['pore.coords'][:,2] += 5.0  # Translate pn2 up 5 units in the Z-direction
        >>> pn.stitch(donor=pn2,pores_1=pn.pores('top'),pores_2=pn2.pores('bottom'),len_max=1.0)
        >>> [pn.Np, pn.Nt]
        [250, 625]
    
        '''
        # Ensure Networks have no associated objects yet
        if (len(network._simulation()) > 1) or (len(donor._simulation()) > 1):
            raise Exception('Cannot stitch a Network with active sibling objects')
        # Get the initial number of pores and throats
        N_init = {}
        N_init['pore']  = network.Np
        N_init['throat'] = network.Nt
        if method == 'delaunay':
            pass
        if method == 'nearest':
            P1 = P_network
            P2 = P_donor + N_init['pore'] # Increment pores on donor
            P = _sp.hstack((P1,P2))
            C1 = network['pore.coords'][P_network]
            C2 = donor['pore.coords'][P_donor]
            D = _sp.spatial.distance.cdist(C1,C2)
            [P1_ind,P2_ind] = _sp.where(D<=len_max)
            conns = _sp.vstack((P1[P1_ind],P2[P2_ind])).T
    
        #Enter donor's pores into the Network
        self.extend(network=network,pore_coords=donor['pore.coords'])
    
        #Enter donor's throats into the Network
        self.extend(network=network,throat_conns=donor['throat.conns']+N_init['pore'])
    
        #Trim throats that are longer then given len_max
        C1 = network['pore.coords'][conns[:,0]]
        C2 = network['pore.coords'][conns[:,1]]
        L = _sp.sum((C1 - C2)**2,axis=1)**0.5
        conns = conns[L<=len_max]
    
        #Add donor labels to recipient network
        if label_suffix != None:
            if label_suffix != '':
                label_suffix = '_'+label_suffix
            for label in donor.labels():
                element = label.split('.')[0]
                locations = _sp.where(network._get_indices(element)>=N_init[element])[0]
                try:
                    network[label+label_suffix]
                except:
                    network[label+label_suffix] = False
                network[label+label_suffix][locations] = donor[label]
    
        #Add the new stitch throats to the Network
        self.extend(network=network,throat_conns=conns,labels='stitched')
        
        # Remove donor from Controller, if present
        # This check allows for the reuse of a donor Network multiple times
        if donor in _ctrl.values():
            _ctrl.purge_object(donor)