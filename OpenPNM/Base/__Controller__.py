'''
###############################################################################
Controller:  Overall simulation controller class
###############################################################################
'''
import pickle as _pickle
import OpenPNM.Base

class Controller(dict):
    r"""

    """
    def __init__(self,obj=None):
        r'''
        '''
        if obj is not None:
            # Update object's simulation attribute
            obj.simulation = self
            # Add object to simulation dict
            self.update({obj.name: obj})

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

    def expand(self):
        r'''
        '''
        obj = list(self.items())[0][1]
        mro = [item.__name__ for item in obj.__class__.__mro__]
        if 'GenericNetwork' not in mro:
            net = obj._net
        else:
            net = obj
        for item in net._geometries:
            item.simulation = self
        for item in net._phases:
            item.simulation = self
        for item in net._physics:
            item.simulation = self

    def clone_object(self,obj):
        r'''
        Clone an OpenPNM Object

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
        '''
        filename = filename.split('.')[0]
        obj = _pickle.load(open(filename+'.pno','rb'))
        obj._sim = self

    def save(self,filename=''):
        r'''
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
        '''
        filename = filename.split('.')[0]
        sim = _pickle.load(open(filename+'.pnm','rb'))
        self.update(sim)
        for item in self.keys():
            self[item]._sim = self

    def export_to(self,filename,fileformat='VTK'):
        r'''
        '''
        pass

    def import_from(self,filename):
        r'''
        '''
        pass


if __name__ == '__main__':
    sim = Controller()












