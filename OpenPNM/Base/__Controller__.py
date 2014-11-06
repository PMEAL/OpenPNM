'''
###############################################################################
Controller:  Overall simulation controller class
###############################################################################
'''

class Controller(dict):
    r"""

    """

    def __init__(self):
        self._add_logger()

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
        '''
        name = obj.name
        # Remove label arrays from all other objects
        for item in self.keys():
            self[item].pop('pore.'+name,None)
            self[item].pop('throat.'+name,None)
        # Set object's simulation to an empty dict
        self[name]._sim = {}
        # Remove object from simulation dict
        del self[name]

    def add(self,obj):
        r'''
        '''
        # Update object's simulation attribute
        obj.simulation = self
        # Add object to simulation dict
        self.update({obj.name: obj})

    def save(self,filename):
        r'''
        '''
        import pickle as _pickle
        filename = filename.split('.')[0]

        #Setup nested dictionary to store simulation
        sim = {}
        #Enter each object's data, object tree and models into dictionary
        for item in self.keys():
            obj = self[item]
            sim[item] = {}
            sim[item]['data'] = obj.copy()
            sim[item]['associations'] = {'Geometries' : obj.geometries(),
                                         'Phases'     : obj.phases(),
                                         'Physics'    : obj.physics()}
            sim[item]['info'] = {'mro'   : [item.__name__ for item in obj.__class__.__mro__],
                                 'class' : obj.__class__}
            # Capture all of the objects attributes, other than OpenPNM reserved ones
            # Note: This could break if people store funky stuff in attributes
            excl_list = ['_physics','_geometries','_phases','_net','_sim', '_name','_logger','_models']
#            a = {key : obj.__dict__[key] for key in obj.__dict__.keys() if key not in excl_list}
#            sim[item]['attrs'] = a
            sim[item]['models'] = {}
            for prop in list(obj._models.keys()):
                sim[item]['models'][prop] = self._save_model(obj,prop)
        #Save nested dictionary pickle
        _pickle.dump(sim,open(filename+'.pnm','wb'))

    def _save_model(self,obj,item):
        r'''
        '''
        #Retrieve info from model
        f = obj._models[item].func
        a = obj._models[item].keywords
        #Store path to model, name of model and argument key:value pairs in a dict
        model = {}
        model['func'] = f
        model['args'] = {}
        for item in a:
            #remove default arguments used in all add_model calls
            if item not in ['physics','network','phase','geometry','mixture']:
                model['args'][item] = a[item]
        return model

    def load(self,filename):
        sim = _pickle.load(open('test.pnm.','rb'))

    @classmethod
    def _add_logger(cls,**kwargs):
        import logging as _logging
        _logging.basicConfig(level=_logging.ERROR,
                             format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                             datefmt='%m-%d %H:%M',
                             )
        if 'loggername' in kwargs.keys():
            cls._logger = _logging.getLogger(kwargs['loggername'])
        else:
            cls._logger = _logging.getLogger(cls.__class__.__name__)
        if 'loglevel' in kwargs.keys():
            cls._loglevel = kwargs['loglevel']


if __name__ == '__main__':
    sim = Controller()












