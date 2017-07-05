# -*- coding: utf-8 -*-
"""
===============================================================================
Utilities.sandbox: A home for functions and classes the we want to audition
===============================================================================

All classes and function in here are highly experimental, subject to change,
and maybe even removal.

"""
import OpenPNM
from OpenPNM.Base import Core
from OpenPNM.Base import logging as _logging
logger = _logging.getLogger(__name__)

class Interface(Core):
    def __init__(self, phase_1, phase_2, **kwargs):
        r"""
        """

        super().__init__(**kwargs)
        if phase_1._net is not phase_2._net:
            raise Exception('Phases are associated with different networks')
        self.update({'pore.all': phase_1['pore.all'].copy()})
        self.update({'throat.all': phase_1['throat.all'].copy()})
        self.phases.update({phase_1.name: phase_1, phase_2.name: phase_2})
        self.network.update({phase_1._net.name: phase_1._net})

        # Add some pore scale models to force interface to have the same
        # temperature and pressure as the phases it separates
        self.models.add(propname='pore.temperature',
                        phase=phase_1,
                        model=OpenPNM.Phases.models.misc.mixture_value)
        self.models.add(propname='pore.pressure',
                        phase=phase_1,
                        model=OpenPNM.Phases.models.misc.mixture_value)

