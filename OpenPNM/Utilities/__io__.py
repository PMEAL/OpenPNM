# -*- coding: utf-8 -*-
from  OpenPNM.Utilities import OpenPNMbase
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
    >>> matfile = OpenPNM.Utilities.ImportMat(filename='example_network',path='D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles')
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
        
class SaveNetwork(OpenPNMbase):
    def __init__(self, **kwargs):
#        self._logger.debug("Write network to file")
        print 'init'
        
    def tocsv(self,net,path='',filename='network'):
        #Write pore info        
        Xp = net.pore_properties['numbering']
        Xp[0] = 'numbering'
        for prop in net.pore_properties.keys():
            if sp.shape(sp.shape(net.pore_properties[prop]))==(1,):
                Xp = sp.vstack((Xp,net.pore_properties[prop]))
        if path=='':
            path = os.path.abspath('')+'\\LocalFiles\\'
        sp.savetxt(path+'\\'+filename+'_pores'+'.csv',Xp.transpose())
        #Write throat info
#        Xt = net.pore_throat['numbering']        
#        for prop in net.throat_properties.keys():
#            Xt = sp.vstack((Xt,net.throat_properties[prop]))
#        sp.savetxt(path+filename+'.csv',Xp)
        
if __name__ == '__main__':
    filename = 'example_network'
    path = 'D:\\AFCC code\\GitHub projects\\OpenPNM\\LocalFiles'
    mat = ImportMat(filename=filename,path=path,loggername="TestImport",loglevel=10)
