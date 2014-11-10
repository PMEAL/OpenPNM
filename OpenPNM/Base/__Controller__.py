'''
###############################################################################
Controller:  Overall simulation controller class
###############################################################################
'''
import pickle as _pickle

class Controller(dict):
    r"""

    """
    # The following __instance class variable and subclassed __new__ method
    # make the Controller class a 'Singleton'.  This way, the _sim attribute
    # of every OpenPNM object is the same, AND if you create a sim on the
    # command line (sim = OpenPNM.Base.Controller()) it will be the same sim!
    __instance = None
    def __new__(cls, *args,**kwargs):
        if Controller.__instance is None:
            Controller.__instance = dict.__new__(cls)
        return Controller.__instance

    def __init__(self):
        r'''
        '''
        pass

    def __str__(self):
        header = ('-'*60)
        print(header)
        print("{a:<25s} {b:<25s}".format(a='Class', b='Object Name'))
        print(header)
        for item in self.keys():
            print("{a:<25s} {b:<25s}".format(a=self[item].__class__.__name__, b=item))
        return ''

    def network(self):
        return self._get_objects(obj_type='GenericNetwork')

    def geometries(self):
        return self._get_objects(obj_type='GenericGeometry')

    def phases(self):
        return self._get_objects(obj_type='GenericPhase')

    def physics(self):
        return self._get_objects(obj_type='GenericPhysics')

    def algorithms(self):
        return self._get_objects(obj_type='GenericAlgorithm')

    def _get_objects(self,obj_type):
        temp = []
        for obj in self.keys():
            mro = [item.__name__ for item in self[obj].__class__.__mro__]
            if obj_type in mro:
                temp.append(self[obj])
        return temp

    def drop(self,obj):
        r'''
        Remove an object from the simulation

        Parameters
        ----------
        obj : OpenPNM Object
            The object to be removed from the simulation.  This method removes
            all traces of the object from everywhere in the simulation,
            including all the object tracking lists and label dictionaries of
            every object.
        '''
        name = obj.name
        for item in self.keys():
            # Remove label arrays from all other objects
            self[item].pop('pore.'+name,None)
            self[item].pop('throat.'+name,None)
            # Remove associations on other objects
            self[item]._geometries[:] = [x for x in self[item]._geometries if x is not obj]
            self[item]._phases[:] = [x for x in self[item]._phases if x is not obj]
            self[item]._physics[:] = [x for x in self[item]._physics if x is not obj]
        # Set object's simulation to an empty dict
        self[name]._sim = {}
        # Remove object from simulation dict
        del self[name]

    def clone_object(self,obj):
        r'''
        Clone an OpenPNM Object, without associating the new object with the
        parent simulation.

        Parameters
        ----------
        obj : OpenPNM Object
            The object to be cloned can be any OpenPNM Object

        Returns
        -------
        A clone of the specified object is returned, but it retains all its links
        to the objects associated with the original object.  The cloned object is
        not associated with the Network.

        Notes
        -----
        This method is intended to create a disposable object, for instance, to
        receive simulation data without overwriting existing data.

        '''
        import pickle
        a = pickle.dumps(obj)
        obj = pickle.loads(a)
        return obj

    def save_object(self,obj,filename=''):
        r'''
        Save a single OpenPNM object to a 'pno' file.

        Parameters
        ----------
        obj : OpenPNM object
            The object to save.  All the associations are removed, so upon
            reloading the object needs to be reconnected manually to a
            simulation.
        filename : string (optional)
            The file name to use when saving.  If no name is given the object
            name is used.
        '''
        if filename == '':
            filename = obj.name
        else:
            filename = filename.split('.')[0]

        obj._sim = {}
        obj._net = []
        obj._geometries = []
        obj._physics = []
        obj._phases = []

        #Save nested dictionary pickle
        _pickle.dump(obj,open(filename+'.pno','wb'))

    def load_object(self,filename):
        r'''
        Load a single object saved as a 'pno' file.

        Parameters
        ----------
        filename : string
            The file name to load.  The file extension must be 'pno'.
        '''
        filename = filename.split('.')[0]
        obj = _pickle.load(open(filename+'.pno','rb'))
        obj._sim = self

    def save(self,filename=''):
        r'''
        Save the entire state of a simulation to a 'pnm' file.

        Parameters
        ----------
        filename : string, optional
            The file name to save as. If none is given the name of the Network
            object is used.
        '''
        if filename == '':
            filename = self.network()[0].name
        else:
            filename = filename.split('.')[0]

        for item in self.keys():
            self[item]._sim = {}

        #Save nested dictionary pickle
        _pickle.dump(self,open(filename+'.pnm','wb'))

    def load(self,filename):
        r'''
        Load an entire simulation from a 'pnm' file.

        Parameters
        ----------
        filename : string
            The file name of the simulation to load.
        '''
        filename = filename.split('.')[0]
        if self != {}:
            print('Warning: Loading data onto non-empty controller object, existing data will be lost')
            self.clear()
        self = _pickle.load(open(filename+'.pnm','rb'))

    def export(self,filename='',fileformat='VTK'):
        r'''
        Export simulation data to the specified file format.

        Parameters
        ----------
        filename : string, optional
            The file name to save as.  If no name is given then the name of
            suppiled object is used.  If no object is given, the name of the
            Network is used.
        obj : OpenPNM object(s), optional
            The object or objects to save.  If no object is given, then the
            entire simulation is saved.
        fileformat : string
            The type of file to create.  Options are:

            1. VTK: Suitable for visualizing in VTK capable software such as Paraview
            2. MAT: Suitable for loading data into Matlab for post-processing

        '''
        if fileformat == 'VTK':
            net = self.network()[0]
            if filename == '':
                filename = net.name
            filename = filename.split('.')[0]
            phases = self.phases()
            import OpenPNM.Utilities.IO as io
            io.VTK.save(filename=filename,network=net,phases=phases)

if __name__ == '__main__':
    sim = Controller()












