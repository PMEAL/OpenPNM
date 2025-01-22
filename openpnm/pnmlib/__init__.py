import numpy as _np
from openpnm.utils import PrintableDict as _pdict
from openpnm.utils import SettingsAttr as _Settings


settings = _Settings()
settings.missing_values = {'bool': False,
                           'int': _np.nan,
                           'float': _np.nan,
                           'object': None}


from . import core
from . import inspect
from . import generators
from . import io
from . import operations
from . import queries
from . import simulations
from . import models
from . import tools
from . import visualization
from . import utils


def info(network):
    r"""
    Prints an overview of the network dictionary
    """
    d = _pdict(network)
    d._key = 'Attribute'
    d._value = 'Description'
    print(d)


def tree(network):
    from openpnm.pnmlib.inspect import tree
    tree(network)
