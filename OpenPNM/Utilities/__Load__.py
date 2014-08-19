import OpenPNM
import scipy as sp
import os, zipfile, shutil

class Load():
    
    @staticmethod
    def simulation(filename):
        r'''
        Load a saved simulation
        '''
        basename = filename.split('.')[0]
        print('Renaming pnm file')
        os.rename(basename+'.pnm',basename+'.zip')
        print('Reading zip archive file')
        zf = zipfile.ZipFile(basename+'.zip')
        #Extract to a temporary directory
        print('Extracting numpy zip files to a temporary directory')
        dirname = 'temp'
        zf.extractall(path=dirname)
        object_list = os.listdir(dirname)
        #Start by initializing the Network
        print('Loading the numpy zip files')
        net = OpenPNM.Network.GenericNetwork(name=dirname)
        obj = sp.load(dirname+'/'+'Network'+'.'+basename+'.npz')
        for item in obj:
            net.update({item:obj[item]})
        #Add Phase, Geometry and Physics objects
        for filename in object_list:
            if filename.split('.')[0] == 'Geometry':
                geom = OpenPNM.Geometry.GenericGeometry(network=net,name=filename.split('.')[1])
                obj = sp.load(dirname+'/'+filename)
                for item in obj:
                    geom.update({item:obj[item]})
            if filename.split('.')[0] == 'Phase':
                phase = OpenPNM.Phases.GenericPhase(network=net,name=filename.split('.')[1])
                obj = sp.load(dirname+'/'+filename)
                for item in obj:
                    phase.update({item:obj[item]})
            if filename.split('.')[0] == 'Physics':
                phase = net.find_object(obj_name=filename.split('.')[3])
                phys = OpenPNM.Physics.GenericPhysics(network=net,phase=phase,name=filename.split('.')[1])
                obj = sp.load(dirname+'/'+filename)
                for item in obj:
                    phys.update({item:obj[item]})
        #Remove the 'temp' directory before exiting
        print('Removing the temporary directory')
        shutil.rmtree(dirname)
        return net
    