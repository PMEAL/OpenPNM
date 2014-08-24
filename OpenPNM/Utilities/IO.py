import OpenPNM
import scipy as _sp
import os as _os

class Save():
    
    @staticmethod
    def PNM(net):
        r'''
        Save the current simulation in it's entirity.  
        
        Parameters
        ----------
        net : OpenPNM Network Object
            The network object of the simulation to be saved.  This will 
            automatically save all Geometry, Phases and Physics objects 
            associated with the Network, but will not save any Algorithms.
        
        '''
        #save simulation as a nested dictionary
        sim = {}
        #Save network
        sim['Network'+'.'+net.name] = net.copy()
        #Save other objects
        for geom in net._geometries:
            sim['Geometry'+'.'+geom.name] = geom.copy()
        for phase in net._phases:
            sim['Phase'+'.'+phase.name] = phase.copy()
            for phys in phase._physics:
                sim['Physics'+'.'+phys.name+'.'+'Phase'+'.'+phase.name] = phys.copy()
        sim['save_info'] = 'just a test'
        _sp.savez_compressed(net.name,**sim)
        _os.rename(net.name+'.npz',net.name+'.pnm')

class Load():
    
    @staticmethod
    def PNM(name):
        r'''
        Load a saved simulation
        
        Parameters
        ----------
        name : string
            The name of the simulation to be read in.
        '''
        #Read in file
        name = name.split('.')[0]
        sim = _sp.load(name+'.pnm')
        #Initializing an empty Network
        net = OpenPNM.Network.GenericNetwork(name=name)
        net.update(sim['Network.'+name].item())
        #Add Phase, Geometry and Physics objects
        for obj in sim.keys():
            if obj.split('.')[0] == 'Geometry':
                geom = OpenPNM.Geometry.GenericGeometry(network=net,name=obj.split('.')[1])
                geom.update(sim[obj].item())
        for obj in sim.keys():
            if obj.split('.')[0] == 'Phase':
                phase = OpenPNM.Phases.GenericPhase(network=net,name=obj.split('.')[1])
                phase.update(sim[obj].item())
        for obj in sim.keys():
            if obj.split('.')[0] == 'Physics':
                phase = net._find_object(obj_name=obj.split('.')[3])
                phys = OpenPNM.Physics.GenericPhysics(network=net,phase=phase,name=obj.split('.')[1])
                phys.update(sim[obj].item())
        return net
    