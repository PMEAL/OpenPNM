from openpnm.core import Base
from openpnm.utils import logging, Docorator
import numpy as np
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.get_sections(base='GenericAlgorithm', sections=['Parameters'])
@docstr.dedent
class GenericAlgorithm(Base):
    r"""
    Generic class to define the foundation of Algorithms.

    Parameters
    ----------
    network : (OpenPNM Network object)
        The network object to which this algorithm will apply.
    name : (string, optional)
        Name of the algorithm
    project : (OpenPNM Project object, optional)
        Either a Network or a Project must be supplied

    Notes
    -----
    This class defines the following methods, which all raise a
    ``NotImplementedError`` and must be defined by the various subclasses

    +---------------------+---------------------------------------------------+
    | Methods             | Description                                       |
    +=====================+===================================================+
    | ``results``         | Generates an array or arrays of data produced by  |
    |                     | the algorithm to be returned to the Phase         |
    +---------------------+---------------------------------------------------+
    | ``setup``           | Collects values to be placed in ``settings``. The |
    |                     | main benefit is defining default values and       |
    |                     | providing documentation on each settings          |
    +---------------------+---------------------------------------------------+
    | ``reset``           | Removes generated data, specified values, and     |
    |                     | any other information lingering on an Algorithm   |
    +---------------------+---------------------------------------------------+

    """

    def __init__(self, network=None, project=None, settings={}, **kwargs):
        self.settings.setdefault('prefix', 'alg')
        self.settings.update(settings)
        if network is not None:
            project = network.project
        super().__init__(project=project, **kwargs)

        # Deal with network or project arguments
        if network is not None:
            if project is not None:
                assert network is project.network
            else:
                project = network.project
        if project:
            self['pore.all'] = np.ones((project.network.Np, ), dtype=bool)
            self['throat.all'] = np.ones((project.network.Nt, ), dtype=bool)

    def results(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")

    def reset(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")

    def setup(self):
        r"""
        """
        raise NotImplementedError("This method must be subclassed")
