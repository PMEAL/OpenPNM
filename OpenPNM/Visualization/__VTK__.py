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
        self._logger.debug("Execute constructor")

    def write(self, net, fluid='none', filename='output.vtp'):
        r"""
        Write Network to a VTK file for visualizing in Paraview

        Parameters
        ----------
        net : OpenPNM Network Object

        filename : string
            Full path to desired file location

        """
        output_path = os.path.join( os.path.expanduser('~'), filename )
        print('     ooooooooooooooooooooooooooooooooooooooooo')
        print('      Writing VTK file:', output_path)
        print('     ooooooooooooooooooooooooooooooooooooooooo')
        self._logger.info("Writing VTK File...please wait")
        self._f = open(output_path,'w')
        self._net=net
        self._fluid=fluid

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
        self._f.write('<Points>\n')
        self._f.write('<DataArray type="Float32" NumberOfComponents="3" format="ascii">\n')
        for i in range(self._net.get_num_pores()):
            self._f.write(str(self._net.pore_properties['coords'][i,0]))
            self._f.write(' ')
            self._f.write(str(self._net.pore_properties['coords'][i,1]))
            self._f.write(' ')
            self._f.write(str(self._net.pore_properties['coords'][i,2]))
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
        pore_keys = list(self._net.pore_properties.keys())
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
                        self._f.write(str(self._net.pore_properties[pore_keys[j]][i]))
                        self._f.write(' ')
                else:
                    for i in range(self._net.get_num_pores()):
                        self._f.write(str(self._net.pore_properties[pore_keys[j]][i][0]))
                        self._f.write(' ')
                self._f.write('\n</DataArray>\n')
        # Now for fluid
        if self._fluid != 'none':
            fluid_name = self._fluid._fluid_recipe['Name']
            pore_keys = list(self._fluid.pore_conditions.keys())
            num_pore_keys = sp.size(pore_keys)
            for j in range(num_pore_keys):
                self._f.write('<DataArray type="Float32" Name="')
                self._f.write(fluid_name+'_'+pore_keys[j])
                self._f.write('" format="ascii">\n')
                size =  np.size(self._fluid.pore_conditions[pore_keys[j]])
                if size == 1:
                    shape =  np.shape(self._fluid.pore_conditions[pore_keys[j]])
                    if np.size(shape) == 0:
                        for i in range(self._net.get_num_pores()):
                            self._f.write(str(np.float(self._fluid.pore_conditions[pore_keys[j]])))
                            self._f.write(' ')
                    else:
                        for i in range(self._net.get_num_pores()):
                            self._f.write(str(np.float(self._fluid.pore_conditions[pore_keys[j]][0])))
                            self._f.write(' ')
                else:
                    shape =  np.shape(self._fluid.pore_conditions[pore_keys[j]])
                    if np.size(shape) == 1:
                        for i in range(self._net.get_num_pores()):
                            self._f.write(str(np.float(self._fluid.pore_conditions[pore_keys[j]][i])))
                            self._f.write(' ')
                    else:
                        for i in range(self._net.get_num_pores()):
                            self._f.write(str(np.float(self._fluid.pore_conditions[pore_keys[j]][i][0])))
                            self._f.write(' ')
                self._f.write('\n</DataArray>\n')
        # Now for fluid.partner
            try:
                fluid_name = self._fluid.partner._fluid_recipe['name']
                pore_keys = self._fluid.partner.pore_conditions.keys()
                num_pore_keys = sp.size(pore_keys)
                for j in range(num_pore_keys):
                    self._f.write('<DataArray type="Float32" Name="')
                    self._f.write(fluid_name+'_'+pore_keys[j])
                    self._f.write('" format="ascii">\n')
                    size =  np.size(self._fluid.partner.pore_conditions[pore_keys[j]])
                    if size == 1:
                        shape =  np.shape(self._fluid.partner.pore_conditions[pore_keys[j]])
                        if np.size(shape) == 0:
                            for i in range(self._net.get_num_pores()):
                                self._f.write(str(np.float(self._fluid.partner.pore_conditions[pore_keys[j]])))
                                self._f.write(' ')
                        else:
                            for i in range(self._net.get_num_pores()):
                                self._f.write(str(np.float(self._fluid.partner.pore_conditions[pore_keys[j]][0])))
                                self._f.write(' ')
                    else:
                        shape =  np.shape(self._fluid.partner.pore_conditions[pore_keys[j]])
                        if np.size(shape) == 1:
                            for i in range(self._net.get_num_pores()):
                                self._f.write(str(np.float(self._fluid.partner.pore_conditions[pore_keys[j]][i])))
                                self._f.write(' ')
                        else:
                            for i in range(self._net.get_num_pores()):
                                self._f.write(str(np.float(self._fluid.partner.pore_conditions[pore_keys[j]][i][0])))
                                self._f.write(' ')
                    self._f.write('\n</DataArray>\n')
            except:
                print('No fluid partner output to VTK')
        self._f.write('</PointData>\n')

    def _write_footer(self):
        self._f.write('</Piece>\n</PolyData>\n</VTKFile>')

