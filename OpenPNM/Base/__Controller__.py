"""
###############################################################################
Controller:  Overall controller class
###############################################################################
"""
import pickle as _pickle
import copy as _copy
import time
import random
import string
import OpenPNM
from OpenPNM.Base import logging
logger = logging.getLogger()


class Controller(dict):
    # The following __instance__ class variable and subclassed __new__ method
    # makes the Controller class a 'Singleton'.  This way, the _ctrl attribute
    # of every OpenPNM object is the same, AND if you create a ctrl on the
    # command line (ctrl = OpenPNM.Base.Controller()) it will be the same ctrl!
    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if Controller.__instance__ is None:
            Controller.__instance__ = dict.__new__(cls)
        return Controller.__instance__

    def __init__(self):
        self.comments = 'Using OpenPNM ' + OpenPNM.__version__

    def __str__(self):
        lines = []
        horizontal_rule = 60 * '-'
        for net in self.networks():
            lines.append(horizontal_rule)
            lines.append('{0:<15} {1:<20} ({2})'.format('Object:',
                                                        'Name',
                                                        'Class'))
            lines.append(horizontal_rule)
            lines.append('{0:<15} {1:<20} ({2})'.format('Network:',
                                                        net.name,
                                                        net.__class__.__name__))
            for geom in net._geometries:
                str = '++ {0:<12} {1:<20} ({2})'
                lines.append(str.format('Geometry: ',
                                        geom.name,
                                        geom.__class__.__name__))
            for phase in net._phases:
                if len(phase._phases) == 0:
                    str = '+ {0:<13} {1:<20} ({2})'
                    lines.append(str.format('Pure Phase: ',
                                            phase.name,
                                            phase.__class__.__name__))
                if len(phase._phases) > 1:
                    str = '+ {0:<13} {1:<20} ({2})'
                    lines.append(str.format('Mixture Phase: ',
                                            phase.name,
                                            phase.__class__.__name__))
                    comps = phase.phases()
                    for compname in comps:
                        str = '++ {0:<12} {1:<20} ({2})'
                        lines.append(str.format('Component Phase: ',
                                                compname,
                                                phase.__class__.__name__))
                for phys in phase._physics:
                    str = '++ {0:<12} {1:<20} ({2})'
                    lines.append(str.format('Physics: ',
                                            phys.name,
                                            phys.__class__.__name__))
        return '\n'.join(lines)

    def _setloglevel(self, level):
        logger.setLevel(level)

    def _getloglevel(self):
        return 'Log level is currently set to: ' + str(logger.level)

    loglevel = property(fget=_getloglevel, fset=_setloglevel)

    def networks(self):
        r"""
        Returns a list of all Network objects
        """
        return self._get_objects(obj_type='GenericNetwork')

    def geometries(self):
        r"""
        Returns a list of all Geometry objects
        """
        return self._get_objects(obj_type='GenericGeometry')

    def phases(self):
        r"""
        Returns a list of all Phase objects
        """
        return self._get_objects(obj_type='GenericPhase')

    def physics(self):
        r"""
        Returns a list of all Physics objects
        """
        return self._get_objects(obj_type='GenericPhysics')

    def algorithms(self):
        r"""
        Returns a list of all Algorithm objects
        """
        return self._get_objects(obj_type='GenericAlgorithm')

    def _get_objects(self, obj_type):
        temp = []
        for obj in self.keys():
            mro = [item.__name__ for item in self[obj].__class__.__mro__]
            if obj_type in mro:
                temp.append(self[obj])
        return temp

    def clear(self):
        r"""
        This is an overloaded version of the standard dict's ``clear`` method.
        This completely clears the Controller object's dict as expected, but
        also removes links to the Controller object in all objects.

        """
        for item in self.keys():
            self[item]._ctrl = {}
        self.__dict__ = {}
        super().clear()

    def update(self, arg):
        r"""
        This is a subclassed version of the standard dict's ``update`` method.
        It can accept a dictionary of OpenPNM Core objects in which case it
        adds the objects to the Controller.  It can also accept an OpenPNM
        Network object in which case it extracts all associated objects and
        adds them to the Controller.  In both cases it adds the Controller to
        all object's ``controller`` attribute.

        Notes
        -----
        The Network (and other Core objects) do not store Algorithms, so this
        update will not add any Algorithm objects to the Controller.  This may
        change.
        """
        if arg.__class__ == dict:
            for item in arg.keys():
                self[item] = arg[item]
                arg[item]._ctrl = self
        else:
            mro = [item.__name__ for item in arg.__class__.__mro__]
            if 'GenericNetwork' in mro:
                net = arg
                self[net.name] = net
                net._ctrl = self
                for item in net._geometries + net._physics + net._phases:
                    self[item.name] = item
                    item._ctrl = self

    def purge_object(self, obj, mode='single'):
        r"""
        Remove an object, including all traces of it in its associated objects

        Parameters
        ----------
        obj : OpenPNM Object
            The object to be removed.  This method removes all traces of the
            object from everywhere, including all the object tracking lists and
            label dictionaries of every object.
        mode : string
            Dicates the type of purge to be performed.  Options are:

            - 'single': Only purges the specified object
            - 'complete': Purges the specified object AND all of its associated
                          objects

        Notes
        -----
        To only remove an object from the Contoller object use the dictionary's
        native ``pop`` method.

        Examples
        --------
        >>> import OpenPNM
        >>> ctrl = OpenPNM.Base.Controller()
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,
        ...                                         pores=pn.Ps,
        ...                                         throats=pn.Ts)

        # Label entries are added to the Network where geom is defined
        >>> 'pore.'+geom.name in pn.keys()
        True
        >>> ctrl.purge_object(geom)

        # geom is removed from Controller object
        >>> geom.name in ctrl.keys()
        False

        # geom's labels are removed from the Network too
        >>> 'pore.' + geom.name in pn.keys()
        False
        """
        if mode == 'complete':
            if obj._net is None:
                net = obj
            else:
                net = obj._net
            for item in net.geometries() + net.phases() + net.physics():
                blank = self.pop(item, None)
            del self[net.name]
        elif mode == 'single':
            name = obj.name
            for item in list(self.keys()):
                # Remove label arrays from all other objects
                self[item].pop('pore.' + name, None)
                self[item].pop('throat.' + name, None)
                # Remove associations on other objects
                self[item]._geometries[:] = \
                    [x for x in self[item]._geometries if x is not obj]
                self[item]._phases[:] = \
                    [x for x in self[item]._phases if x is not obj]
                self[item]._physics[:] = \
                    [x for x in self[item]._physics if x is not obj]
            # Set object's controller attribute to an empty dict
            self[name]._ctrl = {}
            # Remove object from Controller dict
            del self[name]

    def ghost_object(self, obj):
        r"""
        Create a ghost OpenPNM Object containing all the data, methods and
        associations of the original object, but without registering the ghost
        anywhere.   This ghost is intended as a disposable object, for
        instance, to receive data without overwriting existing data.

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
        >>> ctrl = OpenPNM.Base.Controller()
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = ctrl.ghost_object(pn)
        >>> pn is pn2  # A copy of pn is created
        False
        >>> pn2.keys() == pn.keys()  # They have otherwise identical data
        True
        >>> pn2.controller is ctrl # pn2 is not associated with existing Controller
        False

        It can also be used to create ghosts of other object types:

        >>> geom = OpenPNM.Geometry.TestGeometry(network=pn,
        ...                                      pores=pn.Ps,
        ...                                      throats=pn.Ts)
        >>> geo2 = ctrl.ghost_object(geom)
        >>> geom is geo2
        False

        # Ghost has same name as ancestor
        >>> geom.name == geo2.name
        True

        # But they are not the same object
        >>> geo2 is ctrl[geo2.name]
        False

        # The ghost is not registered with the Controller
        >>> geo2.controller
        {}

        # The following comparisons look at some 'behind the scenes' information
        # The ghost and ancestor are assoicated with the same Network
        >>> geo2._net == geom._net
        True

        # But the Network remains aware of the ancestor only
        >>> geo2 in pn._geometries
        False

        """
        obj_new = _copy.copy(obj)
        obj_new.__dict__ = _copy.copy(obj.__dict__)
        obj_new._ctrl = {}
        del self[obj.name]
        self[obj.name] = obj
        return obj_new

    def save_simulation(self, network, filename=''):
        r"""
        Save a single Network simulation to a 'net' file, including all of its
        associated objects, but not Algorithms

        Parameters
        ----------
        network : OpenPNM Network object
            The Network to save
        filename : string, optional
            If no filename is given the name of the Network is used
        """
        if filename == '':
            filename = network.name
        else:
            filename = filename.rstrip('.net')

        network._ctrl = {}
        # Save nested dictionary pickle
        _pickle.dump(network, open(filename + '.net', 'wb'))
        network._ctrl = self

    def load_simulation(self, filename):
        r"""
        Loads a Network simulation fromt the specified 'net' file and adds it
        to the Controller

        Parameters
        ----------
        filename : string
            The name of the file containing the Network simulation to load
        """
        filename = filename.rstrip('.net')
        net = _pickle.load(open(filename + '.net', 'rb'))
        net.controller = self

    def save(self, filename=''):
        r"""
        Save the entire state of the Controller to a 'pnm' file.

        Parameters
        ----------
        filename : string, optional
            The file name to save as. If no filename is provided the current
            date and time is used.

        Examples
        --------
        >>> import OpenPNM
        >>> ctrl = OpenPNM.Base.Controller()
        >>> pn = OpenPNM.Network.TestNet()
        >>> ctrl.save('test.pnm')
        >>> pn.name in ctrl.keys()
        True
        >>> ctrl.clear()
        >>> ctrl.keys()
        dict_keys([])
        >>> ctrl.load('test.pnm')
        >>> pn.name in ctrl.keys()
        True
        """
        if filename == '':
            from datetime import datetime
            i = datetime.now()
            filename = i.strftime('%Y-%m-%d_%H-%M-%S')
        else:
            filename = filename.rstrip('.pnm')

        # Save nested dictionary pickle
        _pickle.dump(self, open(filename + '.pnm', 'wb'))

    def load(self, filename):
        r"""
        Load an entire Controller from a 'pnm' file.

        Parameters
        ----------
        filename : string
            The file name of the Controller to load.

        Notes
        -----
        This calls the ``clear`` method of the Controller object, so it will
        over write the calling objects information AND remove any references
        to the calling object from existing objects.
        """
        filename = filename.strip('.pnm')
        if self != {}:
            print('Warning: Loading data onto non-empty controller object, \
                   existing data will be lost')
            self.clear()
        self = _pickle.load(open(filename + '.pnm', 'rb'))

    def export(self, network=None, filename='', fileformat='VTK'):
        r"""
        Export data to the specified file format.

        Parameters
        ----------
        network : OpenPNM Network Object
            This Network and all of its phases will be written to the specified
            file.  If no Netowrk is given it will check to ensure that only one
            Network exists on the Controller and use that.  If there is more
            than one Network an error is thrown.
        filename : string, optional
            The file name to save as.  If no name is given then the name of
            suppiled object is used.  If no object is given, the name of the
            Network is used.
        fileformat : string
            The type of file to create.  Options are:

            1. VTK: Suitable for visualizing in VTK capable software such as Paraview
            2. MAT: Suitable for loading data into Matlab for post-processing

        """
        if network is None:
            if len(self.networks()) == 1:
                network = self.networks()[0]
            else:
                raise Exception('Multiple Networks found, please specify \
                                which to export')
        import OpenPNM.Utilities.IO as io
        if fileformat == 'VTK':
            phases = network._phases
            io.VTK.save(filename=filename, network=network, phases=phases)
            return
        if fileformat == 'MAT':
            phases = network._phases
            io.MAT.save(filename=filename, network=network, phases=phases)
            return

    def _script(self, filename, mode='read'):
        r"""
        Save or reload the script files used for the modeling

        Parameters
        ----------
        filename : string
            The name of the file to read or write
        mode : string
            Whether to 'archive' the given script file on the object or to
            'retrieve' it from the object and create a new file with it.  The
            default is 'archive'.
        """
        filename = filename.split('.')[0]+'.py'
        if mode == 'archive':
            with open(filename, 'rb') as read_file:
                contents = read_file.read()
            self._script = contents
        if mode == 'retrieve':
            with open(filename, 'wb') as write_file:
                write_file.write(self._script)

    def _set_comments(self, string):
        if hasattr(self, '_comments') is False:
            self._comments = {}
        self._comments[time.strftime('%c')] = string

    def _get_comments(self):
        if hasattr(self, '_comments') is False:
            print('No comments found')
        else:
            for key in self._comments.keys():
                print(key, ': ', self._comments[key])

    comments = property(fget=_get_comments, fset=_set_comments)

    def clone_simulation(self, network, name=None):
        r"""
        Accepts a Network object and creates a complete clone including all
        associated objects.  All objects in the cloned simulation are
        registered with the Controller object and are fully functional.

        Parameters
        ----------
        network : OpenPNM Network Object
            The Network object that is to be cloned.  Because a Network has
            handles to ALL associated objects it acts as the representative
            for the entire simulation.

        name : string
            This string will be appended to the name of all cloned objects.

        Returns
        -------
        A handle to the new Network object, which will include handles to
        clones of all associated objects.

        See Also
        --------
        ghost_object

        Notes
        -----
        One useful application of this method is to create a cloned simulation
        that can be trimmed to a smaller size.  This smaller simulation will
        result in much faster Algorithms calculations.

        Examples
        --------
        >>> import OpenPNM
        >>> ctrl = OpenPNM.Base.Controller()
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = ctrl.clone_simulation(pn, name='cloned')
        >>> pn2 is pn
        False
        """
        if network._parent is not None:
            logger.error('Cannot clone a network that is already a clone')
            return
        if name is None:
            name = ''.join(random.choice(string.ascii_uppercase +
                                         string.ascii_lowercase +
                                         string.digits) for _ in range(5))
        if self._validate_name(network.name + '_' + name) is False:
            logger.error('The provided name is already in use')
            return

        bak = {}
        bak.update(self)
        self.clear()
        net = _copy.deepcopy(network)
        self.update(net)

        # Add supplied name suffix to all cloned objects
        for item in list(self.keys()):
            self[item]._parent = network
            self[item].name = self[item].name + '_' + name

        # Add parent Network numbering to clone
        net['pore.' + network.name] = network.Ps
        net['throat.' + network.name] = network.Ts
        self.update(bak)
        return net

    def _validate_name(self, name):
        valid_name = True
        for item_name in list(self.keys()):
            if name == item_name:
                return False
            for array_name in list(self[item_name].keys()):
                if name == array_name.split('.')[1]:
                    return False
        return valid_name
