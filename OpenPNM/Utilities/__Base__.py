"""
module __OpenPNMbase__: contains OpenPNM base classes
=====================================================

"""
import sys, os, string, random, time
parent_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(1, parent_dir)
import OpenPNM
import scipy.constants
import logging as _logging
import time

# set up logging to file - see previous section for more details
_logging.basicConfig(level=_logging.ERROR,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    )
                    
class Base(object):
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
       DEBUG    10      Detailed information, at diagnostic stage
       INFO     20      Confirmation that things are working as expected.
       WARNING  30      An indication that something unexpected happened.
       ERROR    40      Due to a more serious problem, program might still execute.
       CRITICAL 50      A serious error, might compromise program execution
       ======== =====   =============================================================
       
    """
    _name = None
    def __init__(self,**kwargs):
        super(Base,self).__init__()
        if 'loggername' in kwargs.keys():
            self._logger = _logging.getLogger(kwargs['loggername'])
        else:
            self._logger = _logging.getLogger(self.__class__.__name__)
        if 'loglevel' in kwargs.keys():
            loglevel = kwargs['loglevel']
            self.set_loglevel(loglevel)
        else:
            loglevel = 20
            self.set_loglevel(loglevel)
            
    def __setitem__(self,**kwargs):
        print('Setter not implimented')
              
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
            
    def find_object(self,obj_name='',obj_type=''):
        r'''
        This is a short-cut method.  Given the string name of an 
        OpenPNM Fluid, Geometry, Physics, Algorithm, or Network object 
        this method will return that object
        
        Parameters
        ----------
        obj_name : string
           Name of sought object
        
        Returns
        -------
        OpenPNM Object
        
        
        '''
        if self.__class__.__module__.split('.')[1] == 'Network':
            net = self
        else:
            net = self._net
        
        if obj_name != '':
            objs = []
            if self.name == obj_name:
                return self
            for geom in net._geometries:
                if geom == obj_name:
                    return net._geometries[geom]
                for phys in net._geometries[geom]._physics:
                    if phys == obj_name:
                        return net._geometries[geom]._physics[phys]
                    if net._geometries[geom]._physics[phys]._fluid.name == obj_name:
                        return net._geometries[geom]._physics[phys]._fluid
            return objs # Return empty list if none found
        elif obj_type != '':
            objs = []
            if self.__class__.__module__.split('.')[1] == obj_type:
                return self
            for geom in net._geometries:
                if net._geometries[geom].__class__.__module__.split('.')[1] == obj_type:
                    objs.append(net._geometries[geom])
                for phys in net._geometries[geom]._physics:
                    if net._geometries[geom]._physics[phys].__class__.__module__.split('.')[1] == obj_type:
                        objs.append(net._geometries[geom]._physics[phys])
                    if net._geometries[geom]._physics[phys]._fluid.__class__.__module__.split('.')[1] == obj_type:
                        objs.append(net._geometries[geom]._physics[phys]._fluid)
            return objs

    def _set_name(self,name):
        obj_type = self.__module__.split('.')[1]
        if name == None:
            name = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(20))
            name = obj_type + '_' + name
#        else:
#            self.find_object(obj_name=name)
        self._name = name
    
    def _get_name(self):
        return self._name
        
    name = property(_get_name,_set_name)
    
    def tic(self):
        #Homemade version of matlab tic and toc functions
        global startTime_for_tictoc
        startTime_for_tictoc = time.time()

    def toc(self):
        if 'startTime_for_tictoc' in globals():
            print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
        else:
            print("Toc: start time not set")
            

            
if __name__ == '__main__':
    pass