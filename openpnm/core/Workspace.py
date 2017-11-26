import pickle
import openpnm
import time
from openpnm.core import logging
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

    def _setloglevel(self, level):
        logger.setLevel(level)

    def _getloglevel(self):
        return 'Log level is currently set to: ' + str(logger.level)

    loglevel = property(fget=_getloglevel, fset=_setloglevel)

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
        Loads a simulation from the specified 'net' file and adds it
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
        hr = '―'*80
        s.append(hr)
        for item in self.values():
            s.append(' ' + item.name)
            s.append(item.__str__())
        return '\n'.join(s)
