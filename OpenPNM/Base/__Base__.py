"""
module __OpenPNMbase__: contains OpenPNM base classes
=====================================================

"""
import os, string, random, time, shutil
import OpenPNM
import scipy as sp
import scipy.constants
import logging as _logging


# set up logging to file - see previous section for more details
_logging.basicConfig(level=_logging.ERROR,
                    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                    datefmt='%m-%d %H:%M',
                    )
                    
class Base(dict):
    r"""
    .. class:: `OpenPNM.Base` -- Base class for OpenPNM
    
    
    Base class with a few bells and whistles..
    
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
            loglevel = 30
            self.set_loglevel(loglevel)
        
        #Initialize phase, physics, and geometry tracking lists
        self._phases = []
        self._geometries = []
        self._physics = []
        self._net = []
        
    def __repr__(self):
        return '<%s.%s object at %s>' % (
        self.__class__.__module__,
        self.__class__.__name__,
        hex(id(self)))
              
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
            
    def _find_object(self,obj_name='',obj_type=''):
        r'''
        Find objects associated with a given network model by name or type
        
        Parameters
        ----------
        obj_name : string
           Name of sought object
           
        obj_type : string
            The type of object beign sought.  Options are:
            
            1. 'Network'
            2. 'Geometry'
            3. 'Phases'
            4. 'Physics'
        
        Returns
        -------
        OpenPNM object or list of objects
        
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='geo1')
        >>> temp = pn._find_object(obj_name='geo1')
        >>> temp.name
        'geo1'
        >>> temp = pn._find_object(obj_type='Geometry')
        >>> temp[0].name
        'geo1'
        
        '''
        if self.__module__.split('.')[1] == 'Network':
            net = self
        else:
            net = self._net
        
        if obj_name != '':
            objs = []
            if self.name == obj_name:
                return self
            if net.name == obj_name:
                return net
            for geom in net._geometries:
                if geom.name == obj_name:
                    return geom
            for phase in net._phases:
                if phase.name == obj_name:
                    return phase
            for phys in net._physics:
                if phys.name == obj_name:
                    return phys
            return objs # Return empty list if none found
        elif obj_type != '':
            objs = []
            for geom in net._geometries:
                if geom.__class__.__module__.split('.')[1] == obj_type:
                    objs.append(geom)
            for phys in net._physics:
                if phys.__class__.__module__.split('.')[1] == obj_type:
                    objs.append(phys)
            for phase in net._phases:
                if phase.__class__.__module__.split('.')[1] == obj_type:
                    objs.append(phase)
            return objs
            
    def physics(self,phys_name=[]):
        r'''
        Retrieves Physics associated with the object
        
        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Physics object to retrieve
        Returns
        -------
            If name is NOT provided, then a list of Physics names is returned. 
            If a name or list of names IS provided, then the Physics object(s) 
            with those name(s) is returned.
        '''
        # If arg given as string, convert to list
        if type(phys_name) == str:
            phys_name = [phys_name]
        if phys_name == []:  # If default argument received
            phys = [item.name for item in self._physics]
        else:  # If list of names received
            phys = []
            for item in self._physics:
                if item.name in phys_name:
                    phys.append(item)
        return phys
        
    def phases(self,phase_name=[]):
        r'''
        Retrieves Phases associated with the object
        
        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Phase object(s) to retrieve.  
        Returns
        -------
            If name is NOT provided, then a list of phase names is returned. If
            a name are provided, then a list containing the requested objects 
            is returned.
        '''
        # If arg given as string, convert to list
        if type(phase_name) == str:
            phase_name = [phase_name]
        if phase_name == []:  # If default argument received
            phase = [item.name for item in self._phases]
        else:  # If list of names received
            phase = []
            for item in self._phases:
                if item.name in phase_name:
                    phase.append(item)
        return phase
        
    def geometries(self,geom_name=[]):
        r'''
        Retrieves Geometry object(s) associated with the object
        
        Parameters
        ----------
        name : string or list of strings, optional
            The name(s) of the Geometry object to retrieve.  
        Returns
        -------
            If name is NOT provided, then a list of Geometry names is returned. 
            If a name IS provided, then the Geometry object of that name is 
            returned.
        '''
        # If arg given as string, convert to list
        if type(geom_name) == str:
            geom_name = [geom_name]
        if geom_name == []:  # If default argument received
            geom = [item.name for item in self._geometries]
        else:  # If list of names received
            geom = []
            for item in self._geometries:
                if item.name in geom_name:
                    geom.append(item)
        return geom
        
    def network(self,name=''):
        r'''
        Retrieves the network associated with the object.  If the object is
        a network, then it returns the parent network from which the present
        object derives, or returns an empty list if it has no parents.
        
        Parameters
        ----------
        name : string, optional
            The name of the Geometry object to retrieve.  
            
        Returns
        -------
            If name is NOT provided, then the name of the parent is returned. 
            If a name IS provided, then the parent netowrk object is returned.
            
        Notes
        -----
        This doesn't quite work yet...we have to decide how to treat sub-nets first
        '''
        if name == '':
            try:
                net = self._net.name
            except:
                net = []
        else:
            net = self._net
        return net
            
    def remove_object(self,obj=None,obj_name=''):
        r'''
        Remove specific objects from a model
        
        Parameters
        ----------
        name : string
            The name of the object to delete
            
        Examples
        --------
        >>> pn = OpenPNM.Network.TestNet()
        >>> geom = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='geo')
        >>> geom.name
        'geo'
        >>> pn.delete_object(obj_name='geo')
        >>> pn._find_object(obj_name='geo')
        []
        
        Notes
        -----
        This disassociates the object from the simulation, but does not delete 
        it from memory necessarily.  For instance, the object may still be 
        reachable from the command line.
        
        '''
        if self.__class__.__module__.split('.')[1] == 'Network':
            net = self
        else:
            net = self._net
        if obj_name != '':
            obj = self._find_object(obj_name=obj_name)
        #Get object type, so we know where to delete from
        obj_type = obj.__class__.__module__.split('.')[1]
        if obj_type == 'Geometry':  # Remove geometry from network
            net._geometries.remove(obj)
            net.pop('pore.'+obj.name,None)
            net.pop('throat.'+obj.name,None)
        elif obj_type == 'Phases':
            for phase in net._phases:
                if phase == obj:
                    for physics in phase._physics:
                        physics._phase = None
                        net._physics.remove(physics)
                    net._phases.remove(obj)
                    del phase
        elif obj_type == 'Physics':
            for physics in net._physics:
                if physics == obj:
                    for phase in physics._phase:
                        phase._physics.remove(obj)
                        phase.pop('pore.'+obj.name,None)
                        phase.pop('throat.'+obj.name,None)
                    net._physcis.remove(obj)
                    del phase
                    
    def save(self,filename=''):
        r'''
        '''
        if filename == '':
            filename = self.name
        sp.savez_compressed(filename,**self)
    
    def load(self,filename):
        r'''
        '''
        filename = filename.split('.')[0] + '.npz'
        temp = sp.load(filename)
        return temp
    
    def _set_name(self,name):
        if self._name != None:
            self._logger.error('Renaming objects can have catastrophic consequences')
            return
        if name == None:
            name = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(5))
            name = self.__module__.split('.')[-1].strip('__') + '_' + name
        if self._find_object(obj_name=name) != []:
            self._logger.error('An object with this name already exists')
            return
        self._name = name
    
    def _get_name(self):
        return self._name
        
    name = property(_get_name,_set_name)

if __name__ == '__main__':
    import doctest
    doctest.testmod(verbose=True)



























