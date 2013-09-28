# -*- coding: utf-8 -*-
from  OpenPNM.Base import OpenPNMbase
import scipy as sp
import scipy.io
import numpy as np
import os


         
class ImportMat(OpenPNMbase):
    r"""
    ImportMat - Class for interacting with Matlab files with .mat suffix
    
    Parameters
    ----------
    
    filename : str
        name of file to read from
    path : str
        location of file
    loglevel : int
        Level of the logger (10=Debug, 20=INFO, 30=Warning, 40=Error, 50=Critical)
    
    Examples
    --------
    
    To read in a .mat file with pore network information
    
    >>> import OpenPNM
    >>> matfile = OpenPNM.Base.ImportMat(filename='example_network',path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles')
    >>> pore_volumes = matfile.getvar('pvolumes')
    
  
    """
    def __init__(self,**kwargs):
        r"""
        
        """
        super(ImportMat,self).__init__(**kwargs)
        if 'filename' in kwargs.keys():
            self._filename = kwargs['filename']
        else:
            self._filename = 'example_network'
        if 'path' in kwargs.keys():
            self._path=kwargs['path']
        else:
            self._path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles'
        self._logger.debug("Import from .mat file")
        #Read in the file
        print self._path
        print self._filename
        filepath = os.path.join(self._path,self._filename)
        self._dictionary=scipy.io.loadmat(filepath)
    
    def getvarnames(self):
        r"""
        Returns variable names in the matfile
        """
        keys = self._dictionary.keys()  
        return keys
        
    def getvar(self,varname):
        r"""
        Returns a specific variable from the matfile        
        """
        var = self._dictionary[varname]
        return var

    def getdict(self):
        r"""
        Returns the matfile as a dictionary
        """
        dictionary = self._dictionary
        return dictionary
        
if __name__ == '__main__':
    filename = 'example_network'
    path = 'D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles'
    mat = ImportMat(filename=filename,path=path,loggername="TestImport",loglevel=10)
