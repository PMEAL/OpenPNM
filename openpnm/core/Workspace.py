import pickle
import openpnm
import time
import copy
from openpnm.core import logging
from openpnm.utils import Settings
logger = logging.getLogger()


def singleton(cls):
    r'''
    This is based on PEP318
    '''
    instances = {}

    def getinstance():
        if cls not in instances:
            instances[cls] = cls()
        return instances[cls]
    return getinstance


@singleton
class Workspace(dict):

    def __init__(self):
        super().__init__()
        self.comments = 'Using OpenPNM ' + openpnm.__version__
        self.settings = Settings()

    def _setloglevel(self, level):
        logger.setLevel(level)

    def _getloglevel(self):
        return 'Log level is currently set to: ' + str(logger.level)

    loglevel = property(fget=_getloglevel, fset=_setloglevel)

    def _create_console_handles(self, simulation):
        r"""
        Adds all objects in the given Simulation to the console as variables
        with handle names taken from each object's name.
        """
        import __main__
        for item in simulation:
            __main__.__dict__[item.name] = item

    def save_simulation(self, simulation, filename=''):
        r"""
        Save a single simulation to a 'net' file, including all of its
        associated objects, but not Algorithms

        Parameters
        ----------
        simulation : OpenPNM Network object
            The Network to save

        filename : string, optional
            If no filename is given the name of the Network is used
        """
        if filename == '':
            filename = simulation.name
        else:
            filename = filename.rsplit('.net', 1)[0]

        # Save nested dictionary pickle
        pickle.dump(simulation, open(filename + '.net', 'wb'))

    def load_simulation(self, filename):
        r"""
        Loads a Simulation from the specified 'net' file and adds it
        to the Workspace

        Parameters
        ----------
        filename : string
            The name of the file containing the Network simulation to load
        """
        filename = filename.rsplit('.net', 1)[0]
        sim = pickle.load(open(filename + '.net', 'rb'))
        if sim.name not in self.keys():
            self[sim.name] = sim
        else:
            raise Exception('A simulation with that name is already present')

    def copy_simulation(self, simulation, new_name=None):
        r"""
        Make a copy of an existing Simulation object

        Parameters
        ----------
        simulation : Simulation object
            The Simulation object to be copied

        name : string, optional
            A name for the new copy of the Simulation.  If not supplied, then
            one will be generated (e.g. 'sim_02')

        Returns
        -------
        The new Simulation object

        """
        new_sim = copy.deepcopy(simulation)
        new_sim._name = hex(id(new_sim))  # Assign temporary name
        new_sim.name = new_name
        return new_sim

    def new_simulation(self, name=None):
        r"""
        Creates a new, empty simulation object

        Parameters
        ----------
        name : string (optional)
            The unique name to give to the Simulation.  If none is given, one
            will be automatically generated (e.g. 'sim_01`)

        Returns
        -------
        An empty Simulation object, suitable for passing into a Network
        generator

        """
        sim = openpnm.core.Simulation(name=name)
        return sim

    def import_data(self, filename):
        r"""
        """
        fileformat = filename.split('.')[-1]
        if fileformat.lower() == 'yaml':
            import networkx as nx
            obj = nx.read_yaml(filename)
            openpnm.io.NetworkX.load(obj)

    def export_data(self, simulation, filename=None, fileformat='vtp'):
        r"""
        """
        if filename is None:
            filename = simulation.name + '_' + time.strftime('%Y%b%d_%H%M%p')
        if fileformat.lower() == 'vtp':
            openpnm.io.VTK.save(simulation, filename=filename)
        if fileformat.lower() == 'csv':
            openpnm.io.CSV.save(simulation, filename=filename)
        if fileformat.lower() == 'yaml':
            import networkx as nx
            obj = openpnm.io.NetworkX.save(simulation)
            nx.write_yaml(obj, filename)
        if fileformat.lower() == 'networkx':
            obj = openpnm.io.NetworkX.save(simulation)
            return obj
        if fileformat.lower() == 'mat':
            openpnm.io.MAT.save(simulation, filename=filename)

    def _gen_name(self):
        r"""
        Generates a valid name for simulations
        """
        n = [0]
        for item in self.keys():
            if item.startswith('sim_'):
                n.append(int(item.split('sim_')[1]))
        name = 'sim_'+str(max(n)+1).zfill(2)
        return name

    def _set_comments(self, string):
        if hasattr(self, '_comments') is False:
            self._comments = {}
        self._comments[time.strftime('%c')] = string

    def _get_comments(self):
        if hasattr(self, '_comments') is False:
            self._comments = {}
        for key in list(self._comments.keys()):
            print(key, ': ', self._comments[key])

    comments = property(fget=_get_comments, fset=_set_comments)

    def __str__(self):
        s = []
        hr = '―'*78
        s.append(hr)
        for item in self.values():
            s.append(' ' + item.name)
            s.append(item.__str__())
        return '\n'.join(s)
