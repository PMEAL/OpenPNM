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
        sim['data'] = {}
        sim['tree'] = {}
        sim['mods'] = {}
        
        #Collect all objects into a single list
        all_objs = [net]
        all_objs.extend(net._geometries)
        all_objs.extend(net._phases)
        all_objs.extend(net._physics)
        
        #Enter each object's data, object tree and models into dictionary
        for obj in all_objs:
            module = obj.__module__.split('.')[1]
            sim['data'][module+'.'+obj.name] = obj.copy()
            sim['tree'][module+'.'+obj.name] = {'Geometries' : obj.geometries(),
                                                'Phases'     : obj.phases(),
                                                'Physics'    : obj.physics()}
            sim['mods'][module+'.'+obj.name] = {}
            for prop in list(obj._models.keys()):
                sim['mods'][module+'.'+obj.name][prop] = Save._save_model(obj,prop)
        #Save nested dictionary as a Numpy zip file
        _sp.savez_compressed(net.name,**sim)
        #Rename the zip extension to pnm for kicks
        _os.rename(net.name+'.npz',net.name+'.pnm')
        
    def _save_model(obj,item):
        r'''
        '''
        #Retrieve info from model
        f = obj._models[item].func
        a = obj._models[item].keywords
        #Store path to model, name of model and argument key:value pairs in a dict
        model = {}
        model['propname'] = item
        model['path'] = f.__module__
        model['name'] = f.__name__
        model['args'] = {}
        for item in a:
            #remove default arguments used in all add_model calls
            if item not in ['physics','network','phase','geometry']:
                model['args'][item] = a[item]
        return model

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
        temp = _sp.load(name+'.pnm')
        sim = {}
        sim['data'] = temp['data'].item()
        sim['tree'] = temp['tree'].item()
        sim['mods'] = temp['mods'].item()
        temp.close()

        for obj in sim['data'].keys():  # Network object
            if obj.split('.')[0] == 'Network':
                net = OpenPNM.Network.GenericNetwork(name=obj.split('.')[1])
                net.update(sim['data'][obj])
        
        for obj in sim['data'].keys():  # Geometry objects
            if obj.split('.')[0] == 'Geometry':
                Ps = net.pores(obj.split('.')[1])
                Ts = net.throats(obj.split('.')[1])
                geom = OpenPNM.Geometry.GenericGeometry(network=net,pores=Ps,throats=Ts,name=obj.split('.')[1])
                geom.update(sim['data'][obj])
                for model in sim['mods'][obj].keys():
                    Load._load_model(geom,sim['mods'][obj][model])
        
        for obj in sim['data'].keys():  # Do Pure phases or independent mixtures first
            if (obj.split('.')[0] == 'Phases') and (sim['tree'][obj]['Phases'] == []):
                phase = OpenPNM.Phases.GenericPhase(network=net,name=obj.split('.')[1])
                phase.update(sim['data'][obj])
                for model in sim['mods'][obj].keys():
                    Load._load_model(phase,sim['mods'][obj][model])
        
        for obj in sim['data'].keys():  # Then do proper mixtures which have subphases
            if (obj.split('.')[0] == 'Phases') and (sim['tree'][obj]['Phases'] != []):
                comps = net.phases(sim['tree'][obj]['Phases'])
                #Instantiate mixture phase with list of components
                phase = OpenPNM.Phases.GenericPhase(network=net,name=obj.split('.')[1],components=comps)
                phase.update(sim['data'][obj])
                for model in sim['mods'][obj].keys():
                    Load._load_model(phase,sim['mods'][obj][model])
        
        for obj in sim['data'].keys():  # Physics objects associated with mixures
            if obj.split('.')[0] == 'Physics':
                phase = net.phases(sim['tree'][obj]['Phases'])[0]  # This will always be only 1 phase
                Ps = phase.pores(obj.split('.')[1])
                Ts = phase.throats(obj.split('.')[1])
                phys = OpenPNM.Physics.GenericPhysics(network=net,phase=phase,pores=Ps,throats=Ts,name=obj.split('.')[1])
                phys.update(sim['data'][obj])
                for model in sim['mods'][obj].keys():
                    Load._load_model(phys,sim['mods'][obj][model])
        
        return net
        
    def _load_model(obj,model):
        r'''
        '''
        #Import model using stored path and name
        mod = eval(model['path']+'.'+model['name'])
        #Apply model to object using info in dict
        obj.add_model(model=mod,propname=model['propname'],**model['args'])



    