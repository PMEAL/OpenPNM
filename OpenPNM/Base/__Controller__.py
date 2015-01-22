'''
###############################################################################
Controller:  Overall simulation controller class
###############################################################################
'''
import pickle as _pickle
import copy as _copy
import time, random, string
import OpenPNM
from OpenPNM.Base import logging
logger = logging.getLogger()

class Controller(dict):
    r"""

    """
    # The following __instance__ class variable and subclassed __new__ method
    # makes the Controller class a 'Singleton'.  This way, the _sim attribute
    # of every OpenPNM object is the same, AND if you create a sim on the
    # command line (sim = OpenPNM.Base.Controller()) it will be the same sim!
    __instance__ = None
    def __new__(cls, *args,**kwargs):
        if Controller.__instance__ is None:
            Controller.__instance__ = dict.__new__(cls)
        return Controller.__instance__

    def __init__(self):
        r'''
        '''
        self.comments = 'Using OpenPNM ' + OpenPNM.__version__

    def __str__(self):
        header = ('-'*60)
        print(header)
        print('Networks')
        print(header)
        print("{a:<25s} {b:<25s}".format(a='Class', b='Object Name'))
        print(header)
        for item in self.network():
            print("{a:<25s} {b:<25s}".format(a=item.__class__.__name__, b=item.name))
        print(header)
        print('Geometries')
        print(header)
        print("{a:<25s} {b:<25s}".format(a='Class', b='Object Name'))
        print(header)
        for item in self.geometries():
            print("{a:<25s} {b:<25s}".format(a=item.__class__.__name__, b=item.name))
        print(header)
        print('Phases')
        print(header)
        print("{a:<25s} {b:<25s}".format(a='Class', b='Object Name'))
        print(header)
        for item in self.phases():
            print("{a:<25s} {b:<25s}".format(a=item.__class__.__name__, b=item.name))
        print(header)
        print('Physics')
        print(header)
        print("{a:<25s} {b:<25s}".format(a='Class', b='Object Name'))
        print(header)
        for item in self.physics():
            print("{a:<25s} {b:<25s}".format(a=item.__class__.__name__, b=item.name))
        print(header)
        print('Algorithms')
        print(header)
        print("{a:<25s} {b:<25s}".format(a='Class', b='Object Name'))
        print(header)
        for item in self.algorithms():
            print("{a:<25s} {b:<25s}".format(a=item.__class__.__name__, b=item.name))
        print(header)
        return ''
        
    def _setloglevel(self,level):
        logger.setLevel(level)
        print('Log level has been changed to -->',logger.level)
    
    def _getloglevel(self):
        print('Log level is currently set to -->',logger.level)
        
    loglevel = property(fget=_getloglevel,fset=_setloglevel)

    def show_tree(self):
        r'''
        Prints a heirarchical list of simulation object associations
        '''
        for net in self.network():
            header = ('-'*60)
            print(header)
            print('Network: '+net.name)
            for geom in net.geometries():
                print('+ '+'Geometry: '+geom)
            for phase in net.phases():
                if len(net.phases(phase)[0].phases())==0:
                    print('+ '+'Pure Phase: '+phase)
                if len(net.phases(phase)[0].phases())>1:
                    print('+ '+'Mixture Phase: '+phase)
                    comps = net.phases(phase).phases()
                    for compname in comps:
                        print('+++ '+'Component Phase: '+compname)
                for phys in net.phases(phase)[0].physics():
                    print('++ '+'Physics: '+phys)

    def network(self):
        r'''
        Returns a list of all Network objects in the simulation.
        '''
        return self._get_objects(obj_type='GenericNetwork')

    def geometries(self):
        r'''
        Returns a list of all Geometry objects in the simulation.
        '''
        return self._get_objects(obj_type='GenericGeometry')

    def phases(self):
        r'''
        Returns a list of all Phase objects in the simulation.
        '''
        return self._get_objects(obj_type='GenericPhase')

    def physics(self):
        r'''
        Returns a list of all Physics objects in the simulation.
        '''
        return self._get_objects(obj_type='GenericPhysics')

    def algorithms(self):
        r'''
        Returns a list of all Algorithm objects in the simulation.
        '''
        return self._get_objects(obj_type='GenericAlgorithm')

    def _get_objects(self,obj_type):
        temp = []
        for obj in self.keys():
            mro = [item.__name__ for item in self[obj].__class__.__mro__]
            if obj_type in mro:
                temp.append(self[obj])
        return temp

    def clear(self):
        r'''
        This is an overloaded version of the standard dict's ``clear`` method.
        This completely clears the Controller object's dict as expected, but
        also removes links to the Controller object in all simulation objects.

        '''
        for item in self.keys():
            self[item]._sim = {}
        self.__dict__ = {}
        super(Controller,self).clear()
        
    def update(self,arg):
        r'''
        This is a subclassed version of the standard dict's ``update`` method.  It
        can accept a dictionary as usual, but also an OpenPNM Network object, 
        from which all the associated objects are extracted and added to the 
        controller.
        
        Notes
        -----
        The Network (and other Core objects) do not store Algorithms, so this
        update will not add any Algorithm objects to the Controller.
        '''
        if arg.__class__ == dict:
            super(Controller,self).update(arg)
        else:
            mro = [item.__name__ for item in arg.__class__.__mro__]
            if 'GenericNetwork' in mro:
                net = arg
                self[net.name] = net
                net._sim = self
                for item in net._geometries + net._physics + net._phases:
                    self[item.name] = item
                    item._sim = self

    def purge_object(self,obj):
        r'''
        Remove an object from the simulation

        Parameters
        ----------
        obj : OpenPNM Object
            The object to be removed from the simulation.  This method removes
            all traces of the object from everywhere in the simulation,
            including all the object tracking lists and label dictionaries of
            every object.

        Notes
        -----
        To only remove an object from the Contoller object, without purging all
        traces from the simulation, use the dictionary's native ``pop`` method.

        Examples
        --------
        >>> import OpenPNM
        >>> sim = OpenPNM.Base.Controller()
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=pn.Ps,throats=pn.Ts)
        >>> 'pore.'+geom.name in pn.keys()  # Label entries are added to the Network where geom is defined
        True
        >>> sim.purge_object(geom)
        >>> geom.name in sim.keys()  # geom is removed from Controller object
        False
        >>> 'pore.'+geom.name in pn.keys()  # geom's labels are removed from the Network too
        False
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

    def ghost_object(self,obj):
        r'''
        Create a ghost OpenPNM Object containing all the data, methods and 
        associations of the original object, but without registering the ghost
        anywhere.   This ghost is intended as a disposable object, for 
        instance, to receive simulation data without overwriting existing data.

        Parameters
        ----------
        obj : OpenPNM Object
            The object to be cloned can be any OpenPNM Object

        Returns
        -------
        A clone of the specified object is returned, but it retains all its links
        to the objects associated with the original object.  The cloned object is
        not associated with the Network.

        Examples
        --------
        >>> import OpenPNM
        >>> sim = OpenPNM.Base.Controller()
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = sim.ghost_object(pn)
        >>> pn is pn2  # A clone of pn is created
        False
        >>> pn2.keys() == pn.keys()
        True
        >>> pn2.simulation is sim # pn2 is not associated with existing Controller
        False

        '''
        obj_new = _copy.deepcopy(obj)
        obj_new.__dict__ = obj.__dict__
        obj_new.simulation = {}
        del self[obj.name]
        self[obj.name] = obj
        return obj_new

    def save_object(self,obj,filename=''):
        r'''
        Save a single OpenPNM object to a 'pno' file.  The main purpose of this
        is to save a copy of the object's data dictionary to a file.  This is
        useful for saving Algorithm objects that contain simulation results.

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
        Load a single object saved as a 'pno' file.  This object will be a 
        standalone entity so many methods that expect information about 
        associations will fail, but the numerical data in the dictionary can 
        be accessed.

        Parameters
        ----------
        filename : string
            The file name to load.  The file extension must be 'pno'.
        '''
        filename = filename.split('.')[0]
        obj = _pickle.load(open(filename+'.pno','rb'))
        obj.simulation = self

    def save(self,filename=''):
        r'''
        Save the entire state of a simulation to a 'pnm' file.

        Parameters
        ----------
        filename : string, optional
            The file name to save as. If none is given the name of the Network
            object is used.

        Examples
        --------
        >>> import OpenPNM
        >>> sim = OpenPNM.Base.Controller()
        >>> pn = OpenPNM.Network.TestNet()
        >>> sim.save('test.pnm')
        >>> pn.name in sim.keys()
        True
        >>> sim.clear()
        >>> sim.keys()
        dict_keys([])
        >>> sim.load('test.pnm')
        >>> pn.name in sim.keys()
        True
        '''
        if filename == '':
            filename = self.network()[0].name
        else:
            filename = filename.split('.')[0]

        #Save nested dictionary pickle
        _pickle.dump(self,open(filename+'.pnm','wb'))

    def load(self,filename):
        r'''
        Load an entire simulation from a 'pnm' file.

        Parameters
        ----------
        filename : string
            The file name of the simulation to load.

        Notes
        -----
        This calls the ``clear`` method of the Controller object, so it will
        over write the calling objects information AND remove any references
        to the calling object from existing simulation objects.
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
        fileformat : string
            The type of file to create.  Options are:

            1. VTK: Suitable for visualizing in VTK capable software such as Paraview
            2. MAT: Suitable for loading data into Matlab for post-processing

        '''
        import OpenPNM.Utilities.IO as io
        if fileformat == 'VTK':
            net = self.network()[0]
            phases = net._phases
            io.VTK.save(filename=filename,network=net,phases=phases)
            return
        if fileformat == 'MAT':
            net = self.network()[0]
            phases = net._phases
            io.MAT.save(filename=filename,network=net,phases=phases)
            return

    def _script(self,filename,mode='read'):
        r'''
        Save or reload the script files used for the simulations

        Parameters
        ----------
        filename : string
            The name of the file to read or write
        mode : string
            Whether to 'archive' the given script file on the object or to
            'retrieve' it from the object and create a new file with it.  The
            default is 'archive'.
        '''
        filename = filename.split('.')[0]+'.py'
        if mode == 'archive':
            with open(filename, "rb") as read_file:
                contents = read_file.read()
            self._script = contents
        if mode == 'retrieve':
            with open(filename, "wb") as write_file:
                write_file.write(self._script)

    def _set_comments(self,string):
        if hasattr(self,'_comments') is False:
            self._comments = {}
        self._comments[time.strftime("%c")] = string

    def _get_comments(self):
        if hasattr(self,'_comments') is False:
            print('No comments found')
        else:
            for key in self._comments.keys():
                print(key, ': ', self._comments[key])

    comments = property(fget=_get_comments,fset=_set_comments)

    def clone_simulation(self,network):
        r'''
        Accepts a Network object and creates a complete clone including all 
        associated objects.  The clone simultion is registered with the
        Controller object
        
        Returns
        -------
        A handle to the new network object, which will include handles to 
        clones of all associated objects.  

        See Also
        --------
        clone_object

        Notes
        -----
        The objects in the returned dictionary can be used for simulations as
        usual, but as they are not associated with a Controller, there is
        limited administrative control over them (i.e. saving and such).

        Examples
        --------
        None yet
        '''
        
        bak = {}
        bak.update(self)
        self.clear()
        net = _copy.deepcopy(network)
        self.update(net)
        rand_str = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(5))
        for item in list(self.keys()):
            temp = self.pop(item)
            new_name = temp.name + '_' + rand_str
            temp.name = new_name
            self[temp.name] = temp
        self.update(bak)
        return net

if __name__ == '__main__':
    sim = Controller()












