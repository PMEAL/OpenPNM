'''
###############################################################################
Base:  Abstract Class
###############################################################################
'''
import string, random, collections
import OpenPNM
import scipy as sp
import scipy.constants

class Base(dict):
    r"""
    The abstract Base class for all OpenPNM objects

    Parameters
    ----------
    name : string
        The name for the object.  This must be unique so no two objects in the
        simulation have the same name.  If no name is provided, and random
        string is appended to the objects module name.

    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)

    loggername : string
        Name of the logger. The default is the name of the class.

    """
    _name = None
    _loglevel = 30

    def __new__(typ, *args, **kwargs):
        obj = dict.__new__(typ, *args, **kwargs)
        #Initialize phase, physics, and geometry tracking lists
        obj._phases = []
        obj._geometries = []
        obj._physics = []
        #Initialize ordered dict for storing property models
        obj._models = collections.OrderedDict()
        return obj

    @classmethod
    def _load(cls,name,data):
        r'''
        '''
        inst = cls.__new__(cls)
        inst.update(data)
        inst.name = name
        return inst

    @classmethod
    def _add_logger(cls,**kwargs):
        import logging as _logging
        # set up logging to file - see previous section for more details
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

    def __init__(self,**kwargs):
        super(Base,self).__init__()
        self._add_logger(**kwargs)

    def __repr__(self):
        return '<%s.%s object at %s>' % (
        self.__class__.__module__,
        self.__class__.__name__,
        hex(id(self)))

    def _set_loglevel(self,level=50):
        if type(level) is str:
            desc = {}
            desc['DEBUG']    = 10
            desc['INFO']     = 20
            desc['WARNING']  = 30
            desc['ERROR']    = 40
            desc['CRITICAL'] = 50
            level = string.ascii_uppercase(level)
            level = desc[level]
        self._loglevel = level
        self._logger.setLevel(level)
        self._logger.debug("Changed log level")

    def _get_loglevel(self):
        level = self._loglevel
        desc = {}
        desc[10] = 'DEBUG: Detailed information for for diagnostics and development'
        desc[20] = 'INFO: Confirmation that things are working as expected'
        desc[30] = 'WARNING: An indication that something unexpected happened'
        desc[40] = 'ERROR: Due to a more serious problem, program might still execute'
        desc[50] = 'CRITICAL: A serious error that might compromise program execution'
        print(desc[level])
        return level

    loglevel = property(fget=_get_loglevel,fset=_set_loglevel)

    def set_loglevel(self,level=50):
        r"""
        Sets the effective log level for this class

        Parameters
        ----------
        level : int
            Level above which messages should be logged

        Examples
        --------
        >>> baseobject = OpenPNM.Base.Base()
        >>> baseobject.set_loglevel(30)

        """
        self.loglevel = level

    def _find_object(self,obj_name='',obj_type=''):
        r'''
        Find objects associated with a given network model by name or type

        Parameters
        ----------
        obj_name : string
           Name of sought object

        obj_type : string
            The type of object beign sought.  Options are:

            1. 'Network'
            2. 'Geometry'
            3. 'Phases'
            4. 'Physics'

        Returns
        -------
        OpenPNM object or list of objects

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='geo1',pores=pn.Ps,throats=pn.Ts)
        >>> temp = pn._find_object(obj_name='geo1')
        >>> temp.name
        'geo1'
        >>> temp = pn._find_object(obj_type='Geometry')
        >>> temp[0].name
        'geo1'

        '''
        mro = [item.__name__ for item in self.__class__.__mro__]
        if 'GenericNetwork' in mro:
            net = self
        else:
            net = self._net

        if obj_name != '':
            objs = []
            if self.name == obj_name:
                return self
            if net.name == obj_name:
                return net
            for geom in net._geometries:
                if geom.name == obj_name:
                    return geom
            for phase in net._phases:
                if phase.name == obj_name:
                    return phase
            for phys in net._physics:
                if phys.name == obj_name:
                    return phys
            return objs # Return empty list if none found
        elif obj_type != '':
            obj_type = 'Generic'+obj_type
            objs = []
            for geom in net._geometries:
                if obj_type in [item.__name__ for item in geom.__class__.__mro__]:
                    objs.append(geom)
            for phys in net._physics:
                if obj_type  in [item.__name__ for item in phys.__class__.__mro__]:
                    objs.append(phys)
            for phase in net._phases:
                if obj_type  in [item.__name__ for item in phase.__class__.__mro__]:
                    objs.append(phase)
            return objs

    def physics(self,phys_name=[]):
        r'''
        Retrieves Physics associated with the object

        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Physics object to retrieve
        Returns
        -------
            If name is NOT provided, then a list of Physics names is returned.
            If a name or list of names IS provided, then the Physics object(s)
            with those name(s) is returned.
        '''
        # If arg given as string, convert to list
        if type(phys_name) == str:
            phys_name = [phys_name]
        if phys_name == []:  # If default argument received
            phys = [item.name for item in self._physics]
        else:  # If list of names received
            phys = []
            for item in self._physics:
                if item.name in phys_name:
                    phys.append(item)
        return phys

    def phases(self,phase_name=[]):
        r'''
        Retrieves Phases associated with the object

        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Phase object(s) to retrieve.
        Returns
        -------
            If name is NOT provided, then a list of phase names is returned. If
            a name are provided, then a list containing the requested objects
            is returned.
        '''
        # If arg given as string, convert to list
        if type(phase_name) == str:
            phase_name = [phase_name]
        if phase_name == []:  # If default argument received
            phase = [item.name for item in self._phases]
        else:  # If list of names received
            phase = []
            for item in self._phases:
                if item.name in phase_name:
                    phase.append(item)
        return phase

    def geometries(self,geom_name=[]):
        r'''
        Retrieves Geometry object(s) associated with the object

        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Geometry object to retrieve.
        Returns
        -------
            If name is NOT provided, then a list of Geometry names is returned.
            If a name IS provided, then the Geometry object of that name is
            returned.
        '''
        # If arg given as string, convert to list
        if type(geom_name) == str:
            geom_name = [geom_name]
        if geom_name == []:  # If default argument received
            geom = [item.name for item in self._geometries]
        else:  # If list of names received
            geom = []
            for item in self._geometries:
                if item.name in geom_name:
                    geom.append(item)
        return geom

    def network(self,name=''):
        r'''
        Retrieves the network associated with the object.  If the object is
        a network, then it returns the parent network from which the present
        object derives, or returns an empty list if it has no parents.

        Parameters
        ----------
        name : string, optional
            The name of the Geometry object to retrieve.

        Returns
        -------
            If name is NOT provided, then the name of the parent is returned.
            If a name IS provided, then the parent netowrk object is returned.

        Notes
        -----
        This doesn't quite work yet...we have to decide how to treat sub-nets first
        '''
        if name == '':
            try:
                net = self._net.name
            except:
                net = []
        else:
            net = self._net
        return net

    def remove_object(self,obj=None,obj_name=''):
        r'''
        Remove specific objects from a model

        Parameters
        ----------
        name : string
            The name of the object to delete

        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='geo',pores=pn.Ps,throats=pn.Ts)
        >>> geom.name
        'geo'
        >>> pn.remove_object(obj_name='geo')
        >>> pn._find_object(obj_name='geo')
        []

        Notes
        -----
        This disassociates the object from the simulation, but does not delete
        it from memory necessarily.  For instance, the object may still be
        reachable from the command line.

        '''
        mro = [item.__name__ for item in self.__class__.__mro__]
        if 'GenericNetwork' in mro:
            net = self
        else:
            net = self._net
        if obj_name != '':
            obj = self._find_object(obj_name=obj_name)

        #Get mro for self
        mro = [item.__name__ for item in obj.__class__.__mro__]
        if 'GenericGeometry' in mro:
            net._geometries = [item for item in net._geometries if item is not obj]
            net.pop('pore.'+obj.name,None)
            net.pop('throat.'+obj.name,None)
        elif 'GenericPhase' in mro:
            for phase in net._phases:
                if phase is obj: #Found correct phase
                    for physics in phase._physics:
                        physics._phases = []
                    net._phases = [item for item in net._phases if item is not obj]
                for component in phase._phases: #Cull from other phases too
                    if component is obj:
                        for physics in component._physics:
                            physics._phases = []
                        phase._phases = [item for item in phase._phases if item is not obj]
        elif 'GenericPhysics' in mro:
            for physics in net._physics:
                if physics is obj:
                    for phase in physics._phases:
                        phase._physics = [item for item in phase._physics if item is not obj]
                        phase.pop('pore.'+obj.name,None)
                        phase.pop('throat.'+obj.name,None)
                    net._physics = [item for item in net._physics if item is not obj]

    def save(self,filename=''):
        r'''

        Parameters
        ----------
        filename : string
            The filename to contain the saved object data in Numpy zip format (npz)

        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
        >>> pn.save('test_pn')

        >>> gn = OpenPNM.Network.GenericNetwork()
        >>> gn.load('test_pn')

        >>> # Remove newly created file
        >>> import os
        >>> os.remove('test_pn.npz')

        '''
        if filename == '':
            filename = self.name
        obj_dict = {}
        obj_dict['data'] = self.copy()
        obj_dict['info'] = {}
        obj_dict['info']['name'] = self.name
        obj_dict['info']['module'] = self.__module__
        sp.savez_compressed(filename,**obj_dict)

    def load(self,filename):
        r'''
        Loads a previously saved object's data onto new, empty Generic object

        Parameters
        ----------
        filename : string
            The file containing the saved object data in Numpy zip format (npz)

        Examples
        --------
        >>> pn = OpenPNM.Network.Cubic(shape=[3,3,3])
        >>> pn.save('test_pn')

        >>> gn = OpenPNM.Network.GenericNetwork()
        >>> gn.load('test_pn')

        >>> # Remove newly created file
        >>> import os
        >>> os.remove('test_pn.npz')

        '''
        if (self.Np == 0) and (self.Nt == 0):
            filename = filename.split('.')[0] + '.npz'
            temp = sp.load(filename)
            data_dict = temp['data'].item()
            info_dict = temp['info'].item()
            self.update(data_dict)
            self._name = info_dict['name']
            temp.close()
        else:
            raise Exception('Cannot load saved data onto an active object')

    def OpenPNM_methods(self):
        r'''
        List the OpenPNM methods on the object
        '''
        header = '-'*80
        a = []
        [a.append(item) for item in self.__dir__() if (item[0] != '_') and (item not in dict.__dir__({}))]
        a = sorted(a)
        print(header)
        print("{:<25s}:  Docstring Blurb".format('Method Name'))
        print(header)
        for item in a:
            doc = self.__getattribute__(item).__doc__
            try:
                doc = doc.split('\n')[1]
                doc = doc.lstrip()
            except:
                doc = '---'
            print("{:<25s}:  {}".format(item, doc))
        print(header)

    def _set_name(self,name):
        if self._name != None:
            self._logger.error('Renaming objects can have catastrophic consequences')
            return
        if name == None:
            name = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(5))
            name = self.__module__.split('.')[-1].strip('__') + '_' + name
        if self._find_object(obj_name=name) != []:
            self._logger.error('An object with this name already exists')
            return
        self._name = name

    def _get_name(self):
        return self._name

    name = property(_get_name,_set_name)

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)



























