"""
module __GenericVisualization__: Base class to visualize networks
==================================================================

.. warning:: The classes of this module should be loaded through the 'Visualization/__init__.py' file.

"""

import OpenPNM
import scipy as sp
import numpy as np
import os

from .__GenericVisualization__ import GenericVisualization

class VTK(GenericVisualization):
    r"""
    writeVTK - Class for writing a VTK file

    Parameters
    ----------
    pn : OpenPNM Network Object
        The network which is to be written to the VTK file

    filename : String, optional
        The desire location and file name of the output file.  If not specified defaults to 'output.vtk'

    Examples
    --------
    Create and store a basic network.

    >>> import OpenPNM as PNM
    >>> net = PNM.Generators.SimpleCubic(divisions = [40,40,20],shape=[0.,4.,0.,4.,0.,2.]).generate()
    >>> PNM.Visualization.VTK(net=net,filename = 'testvtp.vtp')

    .. note::
    n/a

    """

    def __init__(self,**kwargs):
        r"""
        Initialize
        """
        super(VTK,self).__init__(**kwargs)
#        self._logger.debug("Execute constructor")

    def write(self, net, fluids='none', filename='output.vtp'):
        r"""
        Write Network to a VTK file for visualizing in Paraview

        Parameters
        ----------
        net : OpenPNM Network Object

        filename : string
            Full path to desired file location

        """
        output_path = os.path.join( os.path.expanduser('~'), filename )
        self._file_name = filename
        self._f = open(output_path,'w')
        self._net=net
        if type(fluids)!= sp.ndarray and fluids!='none': 
            fluids = sp.array(fluids,ndmin=1)
        self._fluids = fluids
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
        self._f.write(str(self._net.num_pores()))
        self._f.write('" NumberOfVerts="0" NumberOfLines="')
        #text = str(pn.getNumThroats())
        self._f.write(str(self._net.num_throats()))
        self._f.write('" NumberOfStrips="0" NumberOfPolys="0">\n')

    def _write_vtk_points(self):
        self._f.write('<Points>\n')
        self._f.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for i in list(range(self._net.num_pores())):
            self._f.write(str(self._net.get_pore_data(prop='coords')[i,0]))
            self._f.write(' ')
            self._f.write(str(self._net.get_pore_data(prop='coords')[i,1]))
            self._f.write(' ')
            self._f.write(str(self._net.get_pore_data(prop='coords')[i,2]))
            self._f.write('\n')
        self._f.write('\n</DataArray>\n</Points>\n')

    def _write_vtk_connections(self):
        self._f.write('<Lines>\n<DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for i in list(range(self._net.num_throats())):
            self._f.write(str(self._net.get_throat_data(prop='connections')[i,0]))
            self._f.write(' ')
            self._f.write(str(self._net.get_throat_data(prop='connections')[i,1]))
            self._f.write(' ')
        self._f.write('\n</DataArray>\n<DataArray type="Int32" Name="offsets" format="ascii">\n')
        for i in list(range(self._net.num_throats())):
            self._f.write(str((i+1)*2))
            self._f.write(' ')
        self._f.write('\n</DataArray>\n</Lines>\n')

    def _write_point_data(self):
        network = self._net
        pore_amalgamate = network.amalgamate_pore_data(fluids=self._fluids)
        throat_amalgamate = network.amalgamate_throat_data(fluids=self._fluids)
        total_amalgamate = dict(pore_amalgamate,**throat_amalgamate)
        total_keys = list(total_amalgamate.keys())
        num_keys = sp.size(total_keys)
        self._f.write('<PointData Scalars="_data">\n')
        for j in list(range(num_keys)):
            if total_keys[j] !='pore_coords' and total_keys[j] !='throat_connections':
                if total_keys[j] in pore_amalgamate: element='pores'
                else: element='throats'
                self._f.write('<DataArray type="Float32" Name="')
                self._f.write(total_keys[j])
                self._f.write('" format="ascii">\n')
                shape =  np.shape(total_amalgamate[total_keys[j]])
                if np.size(shape) == 1:
                    total_amalgamate[total_keys[j]] = total_amalgamate[total_keys[j]]*sp.ones(getattr(network,'num_'+element)())
                    for i in list(range(getattr(network,'num_'+element)())):
                        self._f.write(str(total_amalgamate[total_keys[j]][i]))
                        self._f.write(' ')
                else:
                    for i in list(range(getattr(network,'num_'+element)())):
                        self._f.write(str(total_amalgamate[total_keys[j]][i][0]))
                        self._f.write(' ')
                self._f.write('\n</DataArray>\n')

        output_path = os.path.join( os.path.expanduser('~'), self._file_name )
        print('     oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo')
        print('      Writing VTK file:', output_path)
        print('     oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo')
        self._f.write('</PointData>\n')

    def _write_footer(self):
        self._f.write('</Piece>\n</PolyData>\n</VTKFile>')

