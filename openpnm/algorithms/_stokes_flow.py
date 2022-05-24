import logging
import numpy as np
from copy import deepcopy
from openpnm.algorithms import ReactiveTransport
from openpnm.utils import Docorator, SettingsAttr
logger = logging.getLogger(__name__)
docstr = Docorator()

__all__ = ['StokesFlow']


class StokesFlowSettings:
    quantity = 'pore.pressure'
    conductance = 'throat.hydraulic_conductance'


class StokesFlow(ReactiveTransport):
    r"""
    A subclass of GenericLinearTransport to simulate viscous flow.
    """

    def __init__(self, settings=None, **kwargs):
        self.settings = SettingsAttr(StokesFlowSettings, settings)
        if 'name' not in kwargs.keys():
            kwargs['name'] = 'stokes_01'
        super().__init__(settings=deepcopy(self.settings), **kwargs)
