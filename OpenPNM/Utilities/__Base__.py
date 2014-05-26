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
    _instances = []
    _name = None
    def __init__(self,**kwargs):
        super(Base,self).__init__()
        self._instances.append(self) #Track all instances derived from this class for kicks
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
            
    def __del__(self):
        for item in self._instances:
            if self is item:
                self._instances.remove(item)
        if self.__module__.split('.')[1] == 'Geometry':
            for item in self._physics:
                del self._physics[item].geometry
        print('deleting')

    def save_object(self):
        r'''
        This method saves the object and all of its associated objects. 
        This method is being implemented, we are sorry for the inconvenience.
        '''
        print('not implemented')
        
    def load_object(self):
        r'''
        This method loads an object and all of its associated objects.
        This method is being implemented, we are sorry for the inconvenience.
        '''
        print('not implemented')
        
    def save_object_tocsv(self,path='', filename='', p_prop='all',t_prop='all'):
        r'''
        '''
        if path=='':    path = os.path.abspath('')+'\\LocalFiles\\'
        if filename=='':    filename = self.name
        if type(p_prop)==str:
            if p_prop=='all':
                p_temp = True
                p_prop = self._pore_data.keys()               
            elif p_prop=='not': p_temp = False
            else:
                p_prop = sp.array(p_prop,ndmin=1)
                p_temp = True
        else:
             p_prop = sp.array(p_prop,ndmin=1)
             p_temp = True
        if type(t_prop)==str:
            if t_prop=='all':
                t_temp = True
                t_prop = self._throat_data.keys()               
            elif t_prop=='not': t_temp = False
            else:
                t_prop = sp.array(t_prop,ndmin=1)
                t_temp = True
        else:
             t_prop = sp.array(t_prop,ndmin=1)
             t_temp = True        
    
        if p_temp :
            for p in p_prop:
                if sp.shape(sp.shape(self.get_pore_data(prop=p)))==(1,):
                    sp.savetxt(path+'\\'+filename+'_pores_'+p+'.csv',self.get_pore_data(prop=p))
        if t_temp:
            for t in t_prop:
                if sp.shape(sp.shape(self.get_throat_data(prop=t)))==(1,):
                    sp.savetxt(path+'\\'+filename+'_throats_'+t+'.csv',self.get_throat_data(prop=t))
        
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
        print('Pore data dictionaries:')
        for item in self._pore_data.keys():
            print('  '+item)
        print('Pore info dictionaries:')
        for item in self._pore_info.keys():
            print('  '+item)
        print('Throat data dictionaries:')
        for item in self._throat_data.keys():
            print('  '+item)
        print('Throat info dictionaries:')
        for item in self._throat_info.keys():
            print('  '+item)
            
    def find_object_by_name(self,name):
        r'''
        This is a short-cut method.  Given the string name of an 
        OpenPNM Fluid, Geometry, Physics, Algorithm, or Network object 
        this method will return that object
        
        Parameters
        ----------
        name : string
            Unique name of desired object
        
        Returns
        -------
        OpenPNM Object
        
        
        '''
        for item in self._instances:
            if item.name == name:
                obj = item
        return obj

    def find_object_by_type(self,obj_type):
        r'''
        
        Parameters
        ----------
        obj_type : string
            The type of object to found found.  
            Options the module names (e.g. Network, Geometry, etc). 
            These can be found from obj.__module__, 
            and extracted with obj.__module.split('.')[1]
            
        Returns
        -------
        A dict containing the objects of the type requested.
            
        '''
        obj = {}
        for item in self._instances:
            if item.__module__.split('.')[1] == obj_type:
                obj.update({item.name : item})
        return obj

    def _set_name(self,name):
        obj_type = self.__module__.split('.')[1]
        temp = self.find_object_by_type(obj_type)
        if name == None:
            name = ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for _ in range(20))
            name = obj_type + '_' + name
        for item in temp.keys():
            if obj_type == 'Geometry':
                if self.name:
                    raise Exception('Cannot rename a Geometry')
                for item in self._net.labels(pores='all'):
                    if item == name:
                        raise Exception('Pore label '+name+' already exists')
                for item in self._net.labels(throats='all'):
                    if item == name:
                        raise Exception('Throat label '+name+' already exists')
            if item == name:
                raise Exception('A '+obj_type+' Object with the supplied name already exists')
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
            
    def print_methods(self):
        a = dir(self)
        for item in a:
            if item.split('_')[0] != '':
                print(item)
                
    def print_inheritance(self):
        a = type.mro(type(self))
        a.reverse()
        b = []
        for item in a:
            b.append(dir(item))
        b[-1] = list(set(b[-1]).union(set(dir(self))))
        
        #Remove inherited methods
        for i in range(1,len(b)):
            temp = list(set(b[i]).difference(set(b[i-1])))
            temp.sort()
            priv = []
            pub = []
            for item in temp:
                if item[0] == '_':
                    priv.append(item)
                else:
                    pub.append(item)
            print('==========================================================')            
            print(a[i])
            print('==========================================================')
            print('Private Methods')
            print('----------------------------------------------------------')
            for item in priv:
                print(item)
            print('----------------------------------------------------------')
            print('Public Methods')
            print('----------------------------------------------------------')
            for item in pub:
                print(item)
            print('')
            
            
            
            
            
if __name__ == '__main__':
    pass