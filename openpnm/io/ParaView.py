import os
import numpy as np
import openpnm as op
import subprocess
from flatdict import FlatDict
from openpnm.io import GenericIO


class ParaView(GenericIO):
    r"""
    Class for exporting and viewing OpenPNM networks in ParaView.
    """

    @classmethod
    def export_data(cls, network, filename):
        r"""
        Exports an OpenPNM network to a paraview state file.

        Parameters
        ----------
        network : GenericNetwork
            The network containing the desired data.
        filename : str
            Path to saved .vtp file.

        Notes
        -----
        Outputs a pvsm file that can be opened in Paraview. The pvsm file will
        be saved with the same name as .vtp file

        """
        try:
            import paraview.simple
        except ModuleNotFoundError:
            msg = (
                "The paraview python bindings must be installed using "
                "conda install -c conda-forge paraview, however this may require"
                " using a virtualenv since conflicts with other packages are common."
                " This is why it is not explicitly included as a dependency in"
                " porespy."
            )
            raise ModuleNotFoundError(msg)
        paraview.simple._DisableFirstRenderCameraReset()
        file = os.path.splitext(filename)[0]
        x, y, z = np.ptp(network.coords, axis=0)
        if sum(op.topotools.dimensionality(network)) == 2:
            zshape = 0
            xshape = y
            yshape = x
        elif sum(op.topotools.dimensionality(network)) == 3:
            zshape = x
            xshape = z
            yshape = y
        maxshape = max(xshape, yshape, zshape)
        shape = np.array([xshape, yshape, zshape])
        # Create a new 'XML PolyData Reader'
        Path = os.getcwd() + "\\" + file + '.vtp'
        water = op.phases.Water(network=network)
        net_vtp = paraview.simple.XMLPolyDataReader(FileName=[Path])
        p = op.io.Dict.to_dict(network, phases=[water], element=['pore'],
                               flatten=False,
                               categorize_by=['data' 'object'])
        p = FlatDict(p, delimiter=' | ')
        t = op.io.Dict.to_dict(network, phases=[water], element=['throat'],
                               flatten=False,
                               categorize_by=['data' 'object'])
        t = FlatDict(t, delimiter=' | ')
        net_vtp.CellArrayStatus = t.keys()
        net_vtp.PointArrayStatus = p.keys()
        # Get active view
        render_view = paraview.simple.GetActiveViewOrCreate('RenderView')
        # Uncomment following to set a specific view size
        # render_view.ViewSize = [1524, 527]
        # Get layout
        _ = paraview.simple.GetLayout()
        # Show data in view
        net_vtp_display = paraview.simple.Show(
            net_vtp, render_view, 'GeometryRepresentation'
        )
        # Trace defaults for the display properties.
        net_vtp_display.Representation = 'Surface'
        net_vtp_display.ColorArrayName = [None, '']
        net_vtp_display.OSPRayScaleArray = [
            f'network | {network.name} | labels | pore.all']
        net_vtp_display.OSPRayScaleFunction = 'PiecewiseFunction'
        net_vtp_display.SelectOrientationVectors = 'None'
        net_vtp_display.ScaleFactor = (maxshape-1)/10
        net_vtp_display.SelectScaleArray = 'None'
        net_vtp_display.GlyphType = 'Arrow'
        net_vtp_display.GlyphTableIndexArray = 'None'
        net_vtp_display.GaussianRadius = (maxshape-1)/200
        net_vtp_display.SetScaleArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        net_vtp_display.ScaleTransferFunction = 'PiecewiseFunction'
        net_vtp_display.OpacityArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        net_vtp_display.OpacityTransferFunction = 'PiecewiseFunction'
        net_vtp_display.DataAxesGrid = 'GridAxesRepresentation'
        net_vtp_display.PolarAxes = 'PolarAxesRepresentation'
        # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        net_vtp_display.ScaleTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        net_vtp_display.OpacityTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Reset view to fit data
        render_view.ResetCamera()
        # Get the material library
        _ = paraview.simple.GetMaterialLibrary()
        # Update the view to ensure updated data information
        render_view.Update()
        # Create a new 'Glyph'
        glyph = paraview.simple.Glyph(Input=net_vtp, GlyphType='Arrow')
        glyph.OrientationArray = ['POINTS', 'No orientation array']
        glyph.ScaleArray = ['POINTS', 'No scale array']
        glyph.ScaleFactor = 1
        glyph.GlyphTransform = 'Transform2'
        # Properties modified on glyph
        glyph.GlyphType = 'Sphere'
        glyph.ScaleArray = [
            'POINTS', 'network | ' + network.name +  ' | properties | pore.diameter']
        # Show data in view
        glyph_display = paraview.simple.Show(
            glyph, render_view, 'GeometryRepresentation')
        # Trace defaults for the display properties.
        glyph_display.Representation = 'Surface'
        glyph_display.ColorArrayName = [None, '']
        glyph_display.OSPRayScaleArray = 'Normals'
        glyph_display.OSPRayScaleFunction = 'PiecewiseFunction'
        glyph_display.SelectOrientationVectors = 'None'
        glyph_display.ScaleFactor = (maxshape - 1) / 10
        glyph_display.SelectScaleArray = 'None'
        glyph_display.GlyphType = 'Arrow'
        glyph_display.GlyphTableIndexArray = 'None'
        glyph_display.GaussianRadius = (maxshape - 1) / 200
        glyph_display.SetScaleArray = ['POINTS', 'Normals']
        glyph_display.ScaleTransferFunction = 'PiecewiseFunction'
        glyph_display.OpacityArray = ['POINTS', 'Normals']
        glyph_display.OpacityTransferFunction = 'PiecewiseFunction'
        glyph_display.DataAxesGrid = 'GridAxesRepresentation'
        glyph_display.PolarAxes = 'PolarAxesRepresentation'
        # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        glyph_display.ScaleTransferFunction.Points = [-0.97, 0, 0.5, 0, 0.97, 1, 0.5, 0]
        # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        glyph_display.OpacityTransferFunction.Points = [-0.97, 0, 0.5, 0, 0.97, 1, 0.5, 0]
        # Update the view to ensure updated data information
        render_view.Update()
        # Set active source
        paraview.simple.SetActiveSource(net_vtp)
        # Create a new 'Shrink'
        shrink1 = paraview.simple.Shrink(Input=net_vtp)
        # Properties modified on shrink1
        shrink1.ShrinkFactor = 1.0
        # Show data in view
        shrink_display = paraview.simple.Show(
            shrink1, render_view, 'UnstructuredGridRepresentation')
        # Trace defaults for the display properties.
        shrink_display.Representation = 'Surface'
        shrink_display.ColorArrayName = [None, '']
        shrink_display.OSPRayScaleArray = [
            'network | ' + network.name +  ' | labels | pore.all']
        shrink_display.OSPRayScaleFunction = 'PiecewiseFunction'
        shrink_display.SelectOrientationVectors = 'None'
        shrink_display.ScaleFactor = (maxshape-1)/10
        shrink_display.SelectScaleArray = 'None'
        shrink_display.GlyphType = 'Arrow'
        shrink_display.GlyphTableIndexArray = 'None'
        shrink_display.GaussianRadius = (maxshape-1)/200
        shrink_display.SetScaleArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        shrink_display.ScaleTransferFunction = 'PiecewiseFunction'
        shrink_display.OpacityArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        shrink_display.OpacityTransferFunction = 'PiecewiseFunction'
        shrink_display.DataAxesGrid = 'GridAxesRepresentation'
        shrink_display.PolarAxes = 'PolarAxesRepresentation'
        shrink_display.ScalarOpacityUnitDistance = 1.0349360947089783
        # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        shrink_display.ScaleTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        shrink_display.OpacityTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Hide data in view
        paraview.simple.Hide(net_vtp, render_view)
        # Update the view to ensure updated data information
        render_view.Update()
        # Create a new 'Cell Data to Point Data'
        cellDatatoPointData1 = paraview.simple.CellDatatoPointData(Input=shrink1)
        cellDatatoPointData1.CellDataArraytoprocess = t
        # Show data in view
        cell_data_to_point_data_display = paraview.simple.Show(
            cellDatatoPointData1, render_view, 'UnstructuredGridRepresentation'
        )
        # Trace defaults for the display properties.
        cell_data_to_point_data_display.Representation = 'Surface'
        cell_data_to_point_data_display.ColorArrayName = [None, '']
        cell_data_to_point_data_display.OSPRayScaleArray = [
            'network | ' + network.name +  ' | labels | pore.all']
        cell_data_to_point_data_display.OSPRayScaleFunction = 'PiecewiseFunction'
        cell_data_to_point_data_display.SelectOrientationVectors = 'None'
        cell_data_to_point_data_display.ScaleFactor = (maxshape-1)/10
        cell_data_to_point_data_display.SelectScaleArray = 'None'
        cell_data_to_point_data_display.GlyphType = 'Arrow'
        cell_data_to_point_data_display.GlyphTableIndexArray = 'None'
        cell_data_to_point_data_display.GaussianRadius = (maxshape-1)/200
        cell_data_to_point_data_display.SetScaleArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        cell_data_to_point_data_display.ScaleTransferFunction = 'PiecewiseFunction'
        cell_data_to_point_data_display.OpacityArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        cell_data_to_point_data_display.OpacityTransferFunction = 'PiecewiseFunction'
        cell_data_to_point_data_display.DataAxesGrid = 'GridAxesRepresentation'
        cell_data_to_point_data_display.PolarAxes = 'PolarAxesRepresentation'
        cell_data_to_point_data_display.ScalarOpacityUnitDistance = 1.0349360947089783
        # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        cell_data_to_point_data_display.ScaleTransferFunction.Points = [
            1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        cell_data_to_point_data_display.OpacityTransferFunction.Points = [
            1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Hide data in view
        paraview.simple.Hide(shrink1, render_view)
        # Update the view to ensure updated data information
        render_view.Update()
        # Set active source
        paraview.simple.SetActiveSource(shrink1)
        # Set active source
        paraview.simple.SetActiveSource(cellDatatoPointData1)
        # Set active source
        paraview.simple.SetActiveSource(shrink1)
        # Create a new 'Extract Surface'
        extractSurface1 = paraview.simple.ExtractSurface(Input=shrink1)
        # Show data in view
        extract_surface_display = paraview.simple.Show(
            extractSurface1, render_view, 'GeometryRepresentation')
        # Trace defaults for the display properties.
        extract_surface_display.Representation = 'Surface'
        extract_surface_display.ColorArrayName = [None, '']
        extract_surface_display.OSPRayScaleArray = [
            'network | ' + network.name +  ' | labels | pore.all']
        extract_surface_display.OSPRayScaleFunction = 'PiecewiseFunction'
        extract_surface_display.SelectOrientationVectors = 'None'
        extract_surface_display.ScaleFactor = (maxshape-1)/10
        extract_surface_display.SelectScaleArray = 'None'
        extract_surface_display.GlyphType = 'Arrow'
        extract_surface_display.GlyphTableIndexArray = 'None'
        extract_surface_display.GaussianRadius = (maxshape-1)/200
        extract_surface_display.SetScaleArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        extract_surface_display.ScaleTransferFunction = 'PiecewiseFunction'
        extract_surface_display.OpacityArray = [
            'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
        extract_surface_display.OpacityTransferFunction = 'PiecewiseFunction'
        extract_surface_display.DataAxesGrid = 'GridAxesRepresentation'
        extract_surface_display.PolarAxes = 'PolarAxesRepresentation'
        # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        extract_surface_display.ScaleTransferFunction.Points = [
            1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        extract_surface_display.OpacityTransferFunction.Points = [
            1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Hide data in view
        paraview.simple.Hide(shrink1, render_view)
        # Update the view to ensure updated data information
        render_view.Update()
        # create a new 'Tube'
        tube = paraview.simple.Tube(Input=extractSurface1)
        tube.Scalars = [
            'POINTS', 'network |' + network.name +  '| labels | pore.all']
        tube.Vectors = [None, '1']
        tube.Radius = 0.04
        # Set active source
        paraview.simple.SetActiveSource(extractSurface1)
        # Destroy tube
        paraview.simple.Delete(tube)
        del tube
        # Set active source
        paraview.simple.SetActiveSource(shrink1)
        # Set active source
        paraview.simple.SetActiveSource(cellDatatoPointData1)
        # Set active source
        paraview.simple.SetActiveSource(extractSurface1)
        # Create a new 'Tube'
        tube =paraview.simple.Tube(Input=extractSurface1)
        tube.Scalars = [
            'POINTS', 'network |' + network.name + ' | labels | pore.all']
        tube.Vectors = [None, '1']
        tube.Radius = 0.04
        # Properties modified on tube
        tube.Vectors = ['POINTS', '1']
        # Show data in view
        tube_display = paraview.simple.Show(tube, render_view,
                                            'GeometryRepresentation')
        # Trace defaults for the display properties.
        tube_display.Representation = 'Surface'
        tube_display.ColorArrayName = [None, '']
        tube_display.OSPRayScaleArray = 'TubeNormals'
        tube_display.OSPRayScaleFunction = 'PiecewiseFunction'
        tube_display.SelectOrientationVectors = 'None'
        tube_display.ScaleFactor = (maxshape)/10
        tube_display.SelectScaleArray = 'None'
        tube_display.GlyphType = 'Arrow'
        tube_display.GlyphTableIndexArray = 'None'
        tube_display.GaussianRadius = (maxshape)/200
        tube_display.SetScaleArray = ['POINTS', 'TubeNormals']
        tube_display.ScaleTransferFunction = 'PiecewiseFunction'
        tube_display.OpacityArray = ['POINTS', 'TubeNormals']
        tube_display.OpacityTransferFunction = 'PiecewiseFunction'
        tube_display.DataAxesGrid = 'GridAxesRepresentation'
        tube_display.PolarAxes = 'PolarAxesRepresentation'
        # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
        tube_display.ScaleTransferFunction.Points = [-1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
        tube_display.OpacityTransferFunction.Points = [-1, 0, 0.5, 0, 1, 1, 0.5, 0]
        # Hide data in view
        paraview.simple.Hide(extractSurface1, render_view)
        # Update the view to ensure updated data information
        render_view.Update()
        # Saving camera placements for all active views
        # Current camera placement for render_view
        render_view.CameraPosition = [
            (xshape + 1) / 2, (yshape + 1) / 2, 4.3 * np.sqrt(np.sum(shape / 2)**2)
        ]
        render_view.CameraFocalPoint = [(xi+1) / 2 for xi in shape]
        render_view.CameraParallelScale = np.sqrt(np.sum(shape / 2)**2)
        paraview.simple.SaveState(f"{file}.pvsm")

    @classmethod
    def open_paraview(cls, filename):
        r"""
        Opens a paraview state file directly in ParaView.

        Parameters
        ----------
        filename : str
            Path to input state file.

        """
        file = os.path.splitext(filename)[0]
        statefile = f"{file}.pvsm"
        paraview_path = "paraview"
        subprocess.Popen([paraview_path, statefile])
