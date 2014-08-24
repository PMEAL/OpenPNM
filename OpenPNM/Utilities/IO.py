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
        print('Make temporary directory to store numpy zip files')
        dirname = net.name
        try:
            _os.mkdir(dirname)
        except:
            raise Exception('A model with that name already exists')
            
        #Save network
        print('Save numpy zip files for each object')
        _sp.savez_compressed(dirname+'/'+'Network'+'.'+net.name+'.npz',**net)
        #Save other objects
        for geom in net._geometries:
            filename = dirname+'/'+'Geometry'+'.'+geom.name+'.npz'
            _sp.savez_compressed(filename,**geom)
        for phase in net._phases:
            filename = dirname+'/'+'Phase'+'.'+phase.name+'.npz'
            _sp.savez_compressed(filename,**phase)
            for phys in phase._physics:
                filename = dirname+'/'+'Physics'+'.'+phys.name+'.'+'Phase'+'.'+phase.name+'.npz'
                _sp.savez_compressed(filename,**phys)
        
#        #Convert directory into a 'pnm' file and remove directory
#        print('Generate a zip archive from the temporary directory')
#        shutil.make_archive(base_name=net.name, format="zip", root_dir='temp')
#        print('Remove temporary directory')
#        shutil.rmtree('temp')
#        print('Rename zip file to a pnm file')
#        os.rename(net.name+'.zip',net.name+'.pnm')

class Load():
    
    @staticmethod
    def PNM(simname):
        r'''
        Load a saved simulation
        
        Parameters
        ----------
        sim_name : string
            The name of the simulation to be read in.
        '''
#        basename = filename.split('.')[0]
#        print('Renaming pnm file')
#        os.rename(basename+'.pnm',basename+'.zip')
#        print('Reading zip archive file')
#        zf = zipfile.ZipFile(basename+'.zip')
#        #Extract to a temporary directory
#        print('Extracting numpy zip files to a temporary directory')
#        dirname = 'temp'
#        zf.extractall(path=dirname)
        dirname = simname
        basename = simname
        object_list = _os.listdir(dirname)
        #Start by initializing the Network
        print('Loading the numpy zip files')
        net = OpenPNM.Network.GenericNetwork(name=dirname)
        obj = _sp.load(dirname+'/'+'Network'+'.'+basename+'.npz')
        for item in obj:
            net.update({item:obj[item]})
        #Add Phase, Geometry and Physics objects
        for filename in object_list:
            if filename.split('.')[0] == 'Geometry':
                geom = OpenPNM.Geometry.GenericGeometry(network=net,name=filename.split('.')[1])
                obj = _sp.load(dirname+'/'+filename)
                for item in obj:
                    geom.update({item:obj[item]})
            if filename.split('.')[0] == 'Phase':
                phase = OpenPNM.Phases.GenericPhase(network=net,name=filename.split('.')[1])
                obj = _sp.load(dirname+'/'+filename)
                for item in obj:
                    phase.update({item:obj[item]})
            if filename.split('.')[0] == 'Physics':
                phase = net._find_object(obj_name=filename.split('.')[3])
                phys = OpenPNM.Physics.GenericPhysics(network=net,phase=phase,name=filename.split('.')[1])
                obj = _sp.load(dirname+'/'+filename)
                for item in obj:
                    phys.update({item:obj[item]})
#        #Remove the 'temp' directory before exiting
#        print('Removing the temporary directory')
#        shutil.rmtree(dirname)
        return net
    