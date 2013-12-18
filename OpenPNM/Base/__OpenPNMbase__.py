#! /usr/bin/env python

# Author: CEF PNM Team
# License: TBD
# Copyright (c) 2012

#from __future__ import print_function

"""
module __OpenPNMbase__: contains OpenPNM base classes
=====================================================

.. warning:: The classes of this module should be loaded through the 'Base.__init__.py' file.

Logging:
--------

This module defines the logging format.

Examples:
---------
>>> import OpenPNM 
>>> tmp=PNM.Base.testinheritance()

"""


import logging as _logging
import scipy.constants
from pprint import pprint

# set up logging to file - see previous section for more details
_logging.basicConfig(level=_logging.ERROR,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    )

class Utilities(object):
    r"""
    .. class:: `OpenPNM.Utilities.OpenPNMbase` -- Base class for OpenPNM
    
    
    Base class with a few bells and whistles. Mainly output on screen, logging
    and pickling. This is the class from which all other classes should inherit.
    
    Parameters
    ----------    
    
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    loggername : string
        Name of the logger. The default is the name of the class.

    Attributes
    ----------
    
    self._logger : _logging.logger
       This class defines a logger for the class and all classes inheriting from it.
       it supports all settings of the standard logger. It still needs some fine 
       tuning.
           
       ======== =====   =============================================================
       Level    Value   When it is used
       ======== =====   =============================================================
       DEBUG    10      Detailed information, at d.iagnostic stage
       INFO     20      Confirmation that things are working as expected.
       WARNING  30      An indication that something unexpected happened.
       ERROR    40      Due to a more serious problem, program might still execute.
       CRITICAL 50      A serious error, might compromise program execution
       ======== =====   =============================================================
       
    """
    def __init__(self,**kwargs):
        super(Utilities,self).__init__()
        if 'loggername' in kwargs.keys():
            self._logger = _logging.getLogger(kwargs['loggername'])
        else:
            self._logger = _logging.getLogger(self.__class__.__name__)
        if 'loglevel' in kwargs.keys():
            loglevel=kwargs['loglevel']
            self.set_loglevel(loglevel)
            
        self.constants = scipy.constants
        
    def set_loglevel(self,level=50):
        r"""
        Sets the effective log level for this class
        
        Parameters
        ----------
        level : int
            Level above which messages should be logged.
        """
        self._logger.setLevel(level)
        self._logger.debug("Changed log level")
        
    def print_dicts(self):
        for item in self.pore_conditions.__iter__():
            print(item,': ',type(self.pore_conditions[item]))
        