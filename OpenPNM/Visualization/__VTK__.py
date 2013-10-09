"""
module __GenericVisualization__: Base class to visualize networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Visualization/__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import os

class VTK(OpenPNM.Utilities.OpenPNMbase):
    r"""
    writeVTK - Class for writing a VTK file
    
    Parameters
    ----------
    
    Examples
    --------
    >>> print 'nothing yet'
    
    .. note:: 
    n/a
    
    """
    
    def __init__(self,**kwargs):
        r"""
        Initialize
        """
        super(VTK,self).__init__(**kwargs)
        self._logger.debug("Execute constructor")
       
    def write(self, net, filename='default.vtp', scaling_factor=1):
        r"""
        Write Network to a VTK file for visualizing in Paraview
        
        Parameters
        ----------
        net : OpenPNM Network Object
        
        filename : string
            Full path to desired file location
        
        scaling_factor : int, optional
            Not sure what this does
        """
        output_path = os.path.join( os.path.expanduser('~'), filename )

        self._logger.info("Writing VTK File...please wait")        
        self._f = open(output_path,'w')
        self._net=net

        print( self._net )

        self._scaling_factor = scaling_factor
        
        self._write_vtk_header()
        self._write_vtk_points()
        self._write_vtk_connections()
        self._write_point_data()
        self._write_footer()
        self._f.close()
        
    def _write_vtk_header(self):
        self._f.write('<?xml version="1.0"?>\n')
        self._f.write('<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
        self._f.write('<PolyData>\n')
        self._f.write('<Piece NumberOfPoints="')
        #text = str(pn.getNumPores())
        self._f.write(str(self._net.get_num_pores()))
        self._f.write('" NumberOfVerts="0" NumberOfLines="')
        #text = str(pn.getNumThroats())
        self._f.write(str(self._net.get_num_throats()))
        self._f.write('" NumberOfStrips="0" NumberOfPolys="0">\n')
    
    def _write_vtk_points(self):
        sf = self._scaling_factor
        self._f.write('<Points>\n')
        self._f.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for i in range(self._net.get_num_pores()):
            self._f.write(str(self._net.pore_properties['coords'][i,0]*sf))
            self._f.write(' ')
            self._f.write(str(self._net.pore_properties['coords'][i,1]*sf))
            self._f.write(' ')
            self._f.write(str(self._net.pore_properties['coords'][i,2]*sf))
            self._f.write('\n')
        self._f.write('\n</DataArray>\n</Points>\n')
    
    def _write_vtk_connections(self):
        self._f.write('<Lines>\n<DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for i in range(self._net.get_num_throats()):
            self._f.write(str(self._net.throat_properties['connections'][i,0]))
            self._f.write(' ')
            self._f.write(str(self._net.throat_properties['connections'][i,1]))
            self._f.write(' ')
        self._f.write('\n</DataArray>\n<DataArray type="Int32" Name="offsets" format="ascii">\n')
        for i in range(self._net.get_num_throats()):
            self._f.write(str((i+1)*2))
            self._f.write(' ')            
        self._f.write('\n</DataArray>\n</Lines>\n')
    
    def _write_point_data(self):
        sf = self._scaling_factor
        pore_keys = self._net.pore_properties.keys()
        num_pore_keys = sp.size(pore_keys)
        self._f.write('<PointData Scalars="pore_data">\n')
        for j in range(num_pore_keys):
            if pore_keys[j] !='coords':
                self._f.write('<DataArray type="Float32" Name="')
                self._f.write(pore_keys[j])
                self._f.write('" format="ascii">\n')
                shape =  np.shape(self._net.pore_properties[pore_keys[j]])
                if np.size(shape) == 1:
                    for i in range(self._net.get_num_pores()):
                        if np.dtype(self._net.pore_properties[pore_keys[j]][0])=='float64':
                            self._f.write(str(self._net.pore_properties[pore_keys[j]][i]*sf))
                        else:
                            self._f.write(str(self._net.pore_properties[pore_keys[j]][i]))
                        self._f.write(' ')
                else:
                    for i in range(self._net.get_num_pores()):
                        self._f.write(str(self._net.pore_properties[pore_keys[j]][i][0]))
                        self._f.write(' ')
                self._f.write('\n</DataArray>\n')
        self._f.write('</PointData>\n')
        
    def _write_footer(self):
        self._f.write('</Piece>\n</PolyData>\n</VTKFile>')
        
        