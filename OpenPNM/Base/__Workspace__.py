"""
###############################################################################
Workspace:  A class for managing the workspace of all objects
###############################################################################
"""
import dill as _pickle
import copy as _copy
import time
import random
import string
import OpenPNM
from OpenPNM.Base import logging
logger = logging.getLogger()


class Workspace(dict):
    # The following __instance__ class variable and subclassed __new__ method
    # makes the Workspace class a 'Singleton'.  This way, any instantiation
    # of a workspace object anywhere in the code will return the same object.
    __instance__ = None

    def __new__(cls, *args, **kwargs):
        if Workspace.__instance__ is None:
            Workspace.__instance__ = dict.__new__(cls)
        return Workspace.__instance__

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
                if geom in self.values():
                    lines.append(str.format('Geometry: ',
                                            geom.name,
                                            geom.__class__.__name__))
                else:
                    lines.append(str.format('ERROR: ',
                                            geom.name,
                                            'Object Not in Workspace'))
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
                    if phys in self.values():
                        lines.append(str.format('Physics: ',
                                                phys.name,
                                                phys.__class__.__name__))
                    else:
                        lines.append(str.format('ERROR: ',
                                                phys.name,
                                                'Object Not in Workspace'))
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
        for obj in list(self.keys()):
            mro = [item.__name__ for item in self[obj].__class__.__mro__]
            if obj_type in mro:
                temp.append(self[obj])
        return temp

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
        >>> mgr = OpenPNM.Base.Workspace()
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.GenericGeometry(network=pn,
        ...                                         pores=pn.Ps,
        ...                                         throats=pn.Ts)

        # Label entries are added to the Network where geom is defined
        >>> 'pore.'+geom.name in pn.keys()
        True
        >>> mgr.purge_object(geom)

        # geom is removed from Workspace object
        >>> geom.name in mgr.keys()
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
                self.pop(item, None)
            self.pop(net.name, None)
        elif mode == 'single':
            name = obj.name
            for item in list(self.keys()):
                # Remove label arrays from all other objects
                self[item].pop('pore.' + name, None)
                self[item].pop('throat.' + name, None)
                # Remove associations on other objects
                self[item].geometries.pop(name, None)
                self[item].physics.pop(name, None)
                self[item].phases.pop(name, None)
            # Remove object from Workspace dict
            self.pop(name, None)

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
        >>> mgr = OpenPNM.Base.Workspace()
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = mgr.ghost_object(pn)
        >>> pn is pn2  # A copy of pn is created
        False
        >>> pn2.keys() == pn.keys()  # They have otherwise identical data
        True
        >>> pn2 in mgr.values() # pn2 is not associated with existing Workspace
        False

        It can also be used to create ghosts of other object types:

        >>> geom = OpenPNM.Geometry.TestGeometry(network=pn,
        ...                                      pores=pn.Ps,
        ...                                      throats=pn.Ts)
        >>> geo2 = mgr.ghost_object(geom)
        >>> geom is geo2
        False

        # Ghost has same name as ancestor
        >>> geom.name == geo2.name
        True

        # But they are not the same object
        >>> geo2 is mgr[geo2.name]
        False

        # The ghost is not registered with the Workspace
        >>> geo2 in mgr.values()
        False

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
        self.update({obj.name: obj})
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
            filename = filename.rsplit('.net', 1)[0]

        # Save nested dictionary pickle
        _pickle.dump(network, open(filename + '.net', 'wb'))

    def load_simulation(self, filename):
        r"""
        Loads a Network simulation fromt the specified 'net' file and adds it
        to the Workspace

        Parameters
        ----------
        filename : string
            The name of the file containing the Network simulation to load
        """
        filename = filename.rsplit('.net', 1)[0]
        net = _pickle.load(open(filename + '.net', 'rb'))
        temp_dict = {}  # Store objects temporarily to ensure no exceptions
        if net.name not in self.keys():
            temp_dict[net.name] = net
        else:
            raise Exception('A simulation with that name is already present')
        for item in net._phases + net._physics + net._geometries:
            if item.name not in self.keys():
                temp_dict[item.name] = item
            else:
                raise Exception('An object with that name is already present')
        # If no exceptions, then transfer objects to self
        for item in temp_dict.values():
            item.workspace = self

    def save_workspace(self, filename=''):
        r"""
        Save the entire state of the Workspace to a 'pnm' file.

        Parameters
        ----------
        filename : string, optional
            The file name to save as. If no filename is provided the current
            date and time is used.

        Examples
        --------

        .. code-block:: python

            import OpenPNM
            mgr = OpenPNM.Base.Workspace()
            mgr.clear()  # Ensure no previous objects are present
            pn = OpenPNM.Network.TestNet()
            mgr.save('test.pnm')
            pn.name in mgr.keys()
            #=> True
            mgr.clear()
            mgr.keys()
            dict_keys([])
            mgr.load('test.pnm')
            pn.name in mgr.keys()
            #=> True

        """
        if filename == '':
            from datetime import datetime
            i = datetime.now()
            filename = i.strftime('%Y-%m-%d_%H-%M-%S')
        else:
            filename = filename.rstrip('.pnm')

        # Save nested dictionary pickle
        _pickle.dump(self, open(filename + '.pnm', 'wb'))

    def save(self, **kwargs):
        r"""
        This method is deprecated, use ``save_workspace`` instead.
        """
        logger.warning("This method is deprecated, use \'save_workspace\'.")
        self.save_workspace(**kwargs)

    def load_workspace(self, filename):
        r"""
        Load an entire Workspace from a 'pnm' file.

        Parameters
        ----------
        filename : string
            The file name of the Workspace to load.

        Notes
        -----
        This calls the ``clear`` method of the Workspace object, so it will
        remove all existing objects in the current workspace.
        """
        filename = filename.rsplit('.pnm', 1)[0]
        if self != {}:
            logger.warn('Loading data onto non-empty workspace object,' +
                        ' existing data will be lost')
            self.clear()

        self = _pickle.load(open(filename+'.pnm', 'rb'))
        for item in self._comments.values():
            if 'Using OpenPNM' in item:
                version = item.lstrip('Using OpenPNM ')
                if version < OpenPNM.__version__:
                    logger.warning('File was created with an earlier version ' +
                                   'OpenPNM: \n' +
                                   '--> File saved with version: ' +
                                   str(version) +
                                   '\n' +
                                   '--> Current version: ' +
                                   str(OpenPNM.__version__))

    def load(self, **kwargs):
        r"""
        This method is deprecated, use ``load_workspace`` instead.
        """
        logger.warning("This method is deprecated, use \'load_workspace\'.")
        self.load_workspace(**kwargs)

    def export(self, network=None, filename='', fileformat='VTK'):
        logger.warning("This method is deprecated, use \'export_data\'.")
        self.export_data(network=network,
                         filename=filename,
                         fileformat=fileformat)

    def export_data(self, network=None, filename='', fileformat='VTK'):
        r"""
        Export data to the specified file format.

        Parameters
        ----------
        network : OpenPNM Network Object
            This Network and all of its phases will be written to the specified
            file.  If no Network is given it will check to ensure that only one
            Network exists on the Workspace and use that.  If there is more
            than one Network an error is thrown.
        filename : string, optional
            The file name to save as.  If no name is given then the name of
            suppiled object is used.  If no object is given, the name of the
            Network is used.
        fileformat : string
            The type of file to create.  Options are:

            **'VTK'**: Suitable for visualizing in VTK capable software such
            as Paraview

            **'MAT'**: Suitable for loading data into Matlab for post-
            processing

            **'CSV'**: Suitable for analyzing data in a spreadsheet program
            such as Excel.  This will save two files, one containing pore data
            and one containing throat data.  These indicators are appended to
            to file names.

        """
        import OpenPNM.Utilities.IO as io

        if network is None:
            if len(self.networks()) == 1:
                network = self.networks()[0]
            else:
                raise Exception('Multiple Networks found, please specify' +
                                'which to export')
        # Generate filename if necessary
        if filename == '':
            filename = network.name

        fileformat = fileformat.lower()
        if fileformat == 'vtk':
            phases = network._phases
            io.VTK.save(filename=filename, network=network, phases=phases)
            return
        elif fileformat == 'mat':
            phases = network._phases
            io.MAT.save(filename=filename, network=network, phases=phases)
            return
        elif fileformat == 'csv':
            phases = network._phases
            io.CSV.save(network=network, filename=filename, phases=phases)
        else:
            raise ValueError(fileformat+' is not a valid format')

    def import_data(self, filename=None):
        r"""
        Import network data stored in an external file format.

        Parameters
        ----------
        filename : string
            The name of the file containing the data.  This should include the
            file extension, which must be one of the following supported
            formats:

            **'csv'** : Comma-separated values as typically used in spreadsheet
            type programs

            **'mat'** : A Matlab \'mat-file\'

            **'vtp'** : A VTK file format used by programs like Paraview

            **'yaml'** : A NetworkX output format

        Notes
        -----
        This is a wrapper or convenience method for the actual IO classes
        located in OpenPNM.Utilities.IO. Refer to the doc strings for those
        classes for information on the actual file format specifications.

        """
        import OpenPNM.Utilities.IO as io
        # Handle normal file types
        ext = filename.split('.')[-1]
        if ext.lower() == 'csv':
            network = io.CSV.load(filename=filename)
        elif ext.lower() == 'mat':
            network = io.MAT.load(filename=filename)
        elif ext.lower() == 'yaml':
            network = io.NetworkX.load(filename=filename)
        elif ext.lower() == 'vtp':
            network = io.VTK.load(filename=filename)
        else:
            raise Exception('Filename does not have suppored extension')
        return network

    def _set_comments(self, string):
        if hasattr(self, '_comments') is False:
            self._comments = {}
        self._comments[time.strftime('%c')] = string

    def _get_comments(self):
        if hasattr(self, '_comments') is False:
            logger.info('No comments found')
        else:
            for key in list(self._comments.keys()):
                logger.info(key, ': ', self._comments[key])

    comments = property(fget=_get_comments, fset=_set_comments)

    def clone_simulation(self, network, name=None):
        r"""
        Accepts a Network object and creates a complete clone including all
        associated objects.  All objects in the cloned simulation are
        registered with the Workspace object and are fully functional.

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
        >>> mgr = OpenPNM.Base.Workspace()
        >>> pn = OpenPNM.Network.TestNet()
        >>> pn2 = mgr.clone_simulation(pn, name='cloned')
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

        net = _copy.deepcopy(network)  # Make clone
        # Add supplied name suffix to all cloned objects
        for item in net._simulation():
            item._parent = network
            item.name = item.name + '_' + name

        # Add parent Network numbering to clone
        net['pore.' + network.name] = network.Ps
        net['throat.' + network.name] = network.Ts
        return net

    def _validate_name(self, name):
        valid_name = True
        for item_name in list(self.keys()):
            # Check object names for conflict
            if name == item_name:
                return False
            # Also check array names on all objects
            for array_name in list(self[item_name].keys()):
                if name == array_name.split('.')[-1]:
                    return False
        return valid_name
