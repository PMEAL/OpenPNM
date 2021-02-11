import os
import numpy as np
import openpnm as op
import subprocess
from flatdict import FlatDict


def export_data(network, filename):
    r"""
    Converts an array to a paraview state file.

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
    net555vtp = paraview.simple.XMLPolyDataReader(FileName=[Path])
    p = op.io.Dict.to_dict(network, phases=[water], element=['pore'],
                           flatten=False,
                           categorize_by=['data' 'object'])
    p = FlatDict(p, delimiter=' | ')
    t = op.io.Dict.to_dict(network, phases=[water], element=['throat'],
                           flatten=False,
                           categorize_by=['data' 'object'])
    t = FlatDict(t, delimiter=' | ')
    net555vtp.CellArrayStatus = t.keys()
    net555vtp.PointArrayStatus = p.keys()
    # Get active view
    renderView1 = paraview.simple.GetActiveViewOrCreate('RenderView')
    # Uncomment following to set a specific view size
    # renderView1.ViewSize = [1524, 527]
    # Get layout
    layout1 = paraview.simple.GetLayout()
    # Show data in view
    net555vtpDisplay = paraview.simple.Show(
        net555vtp, renderView1, 'GeometryRepresentation'
    )
    # Trace defaults for the display properties.
    net555vtpDisplay.Representation = 'Surface'
    net555vtpDisplay.ColorArrayName = [None, '']
    net555vtpDisplay.OSPRayScaleArray = [f'network | {network.name} | labels | pore.all']
    net555vtpDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
    net555vtpDisplay.SelectOrientationVectors = 'None'
    net555vtpDisplay.ScaleFactor = (maxshape-1)/10
    net555vtpDisplay.SelectScaleArray = 'None'
    net555vtpDisplay.GlyphType = 'Arrow'
    net555vtpDisplay.GlyphTableIndexArray = 'None'
    net555vtpDisplay.GaussianRadius = (maxshape-1)/200
    net555vtpDisplay.SetScaleArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    net555vtpDisplay.ScaleTransferFunction = 'PiecewiseFunction'
    net555vtpDisplay.OpacityArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    net555vtpDisplay.OpacityTransferFunction = 'PiecewiseFunction'
    net555vtpDisplay.DataAxesGrid = 'GridAxesRepresentation'
    net555vtpDisplay.PolarAxes = 'PolarAxesRepresentation'
    # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    net555vtpDisplay.ScaleTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    net555vtpDisplay.OpacityTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Reset view to fit data
    renderView1.ResetCamera()
    # Get the material library
    materialLibrary1 = paraview.simple.GetMaterialLibrary()
    # Update the view to ensure updated data information
    renderView1.Update()
    # Create a new 'Glyph'
    glyph1 = paraview.simple.Glyph(Input=net555vtp, GlyphType='Arrow')
    glyph1.OrientationArray = ['POINTS', 'No orientation array']
    glyph1.ScaleArray = ['POINTS', 'No scale array']
    glyph1.ScaleFactor = 1
    glyph1.GlyphTransform = 'Transform2'
    # Properties modified on glyph1
    glyph1.GlyphType = 'Sphere'
    glyph1.ScaleArray = [
        'POINTS', 'network | ' + network.name +  ' | properties | pore.diameter']
    # Show data in view
    glyph1Display = paraview.simple.Show(glyph1, renderView1, 'GeometryRepresentation')
    # Trace defaults for the display properties.
    glyph1Display.Representation = 'Surface'
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.OSPRayScaleArray = 'Normals'
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.SelectOrientationVectors = 'None'
    glyph1Display.ScaleFactor = (maxshape - 1) / 10
    glyph1Display.SelectScaleArray = 'None'
    glyph1Display.GlyphType = 'Arrow'
    glyph1Display.GlyphTableIndexArray = 'None'
    glyph1Display.GaussianRadius = (maxshape - 1) / 200
    glyph1Display.SetScaleArray = ['POINTS', 'Normals']
    glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
    glyph1Display.OpacityArray = ['POINTS', 'Normals']
    glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
    glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
    glyph1Display.PolarAxes = 'PolarAxesRepresentation'
    # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    glyph1Display.ScaleTransferFunction.Points = [-0.97, 0, 0.5, 0, 0.97, 1, 0.5, 0]
    # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    glyph1Display.OpacityTransferFunction.Points = [-0.97, 0, 0.5, 0, 0.97, 1, 0.5, 0]
    # Update the view to ensure updated data information
    renderView1.Update()
    # Set active source
    paraview.simple.SetActiveSource(net555vtp)
    # Create a new 'Shrink'
    shrink1 = paraview.simple.Shrink(Input=net555vtp)
    # Properties modified on shrink1
    shrink1.ShrinkFactor = 1.0
    # Show data in view
    shrink1Display = paraview.simple.Show(
        shrink1, renderView1, 'UnstructuredGridRepresentation'
    )
    # Trace defaults for the display properties.
    shrink1Display.Representation = 'Surface'
    shrink1Display.ColorArrayName = [None, '']
    shrink1Display.OSPRayScaleArray = [
        'network | ' + network.name +  ' | labels | pore.all']
    shrink1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    shrink1Display.SelectOrientationVectors = 'None'
    shrink1Display.ScaleFactor = (maxshape-1)/10
    shrink1Display.SelectScaleArray = 'None'
    shrink1Display.GlyphType = 'Arrow'
    shrink1Display.GlyphTableIndexArray = 'None'
    shrink1Display.GaussianRadius = (maxshape-1)/200
    shrink1Display.SetScaleArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    shrink1Display.ScaleTransferFunction = 'PiecewiseFunction'
    shrink1Display.OpacityArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    shrink1Display.OpacityTransferFunction = 'PiecewiseFunction'
    shrink1Display.DataAxesGrid = 'GridAxesRepresentation'
    shrink1Display.PolarAxes = 'PolarAxesRepresentation'
    shrink1Display.ScalarOpacityUnitDistance = 1.0349360947089783
    # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    shrink1Display.ScaleTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    shrink1Display.OpacityTransferFunction.Points = [1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Hide data in view
    paraview.simple.Hide(net555vtp, renderView1)
    # Update the view to ensure updated data information
    renderView1.Update()
    # Create a new 'Cell Data to Point Data'
    cellDatatoPointData1 = paraview.simple.CellDatatoPointData(Input=shrink1)
    cellDatatoPointData1.CellDataArraytoprocess = t
    # Show data in view
    cellDatatoPointData1Display = paraview.simple.Show(
        cellDatatoPointData1, renderView1, 'UnstructuredGridRepresentation'
    )
    # Trace defaults for the display properties.
    cellDatatoPointData1Display.Representation = 'Surface'
    cellDatatoPointData1Display.ColorArrayName = [None, '']
    cellDatatoPointData1Display.OSPRayScaleArray = [
        'network | ' + network.name +  ' | labels | pore.all']
    cellDatatoPointData1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    cellDatatoPointData1Display.SelectOrientationVectors = 'None'
    cellDatatoPointData1Display.ScaleFactor = (maxshape-1)/10
    cellDatatoPointData1Display.SelectScaleArray = 'None'
    cellDatatoPointData1Display.GlyphType = 'Arrow'
    cellDatatoPointData1Display.GlyphTableIndexArray = 'None'
    cellDatatoPointData1Display.GaussianRadius = (maxshape-1)/200
    cellDatatoPointData1Display.SetScaleArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    cellDatatoPointData1Display.ScaleTransferFunction = 'PiecewiseFunction'
    cellDatatoPointData1Display.OpacityArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    cellDatatoPointData1Display.OpacityTransferFunction = 'PiecewiseFunction'
    cellDatatoPointData1Display.DataAxesGrid = 'GridAxesRepresentation'
    cellDatatoPointData1Display.PolarAxes = 'PolarAxesRepresentation'
    cellDatatoPointData1Display.ScalarOpacityUnitDistance = 1.0349360947089783
    # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    cellDatatoPointData1Display.ScaleTransferFunction.Points = [
        1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    cellDatatoPointData1Display.OpacityTransferFunction.Points = [
        1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Hide data in view
    paraview.simple.Hide(shrink1, renderView1)
    # Update the view to ensure updated data information
    renderView1.Update()
    # Set active source
    paraview.simple.SetActiveSource(shrink1)
    # Set active source
    paraview.simple.SetActiveSource(cellDatatoPointData1)
    # Set active source
    paraview.simple.SetActiveSource(shrink1)
    # Create a new 'Extract Surface'
    extractSurface1 = paraview.simple.ExtractSurface(Input=shrink1)
    # Show data in view
    extractSurface1Display = paraview.simple.Show(extractSurface1, renderView1,
                                                  'GeometryRepresentation')
    # Trace defaults for the display properties.
    extractSurface1Display.Representation = 'Surface'
    extractSurface1Display.ColorArrayName = [None, '']
    extractSurface1Display.OSPRayScaleArray = [
        'network | ' + network.name +  ' | labels | pore.all']
    extractSurface1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    extractSurface1Display.SelectOrientationVectors = 'None'
    extractSurface1Display.ScaleFactor = (maxshape-1)/10
    extractSurface1Display.SelectScaleArray = 'None'
    extractSurface1Display.GlyphType = 'Arrow'
    extractSurface1Display.GlyphTableIndexArray = 'None'
    extractSurface1Display.GaussianRadius = (maxshape-1)/200
    extractSurface1Display.SetScaleArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    extractSurface1Display.ScaleTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.OpacityArray = [
        'POINTS', 'network | ' + network.name +  ' | labels | pore.all']
    extractSurface1Display.OpacityTransferFunction = 'PiecewiseFunction'
    extractSurface1Display.DataAxesGrid = 'GridAxesRepresentation'
    extractSurface1Display.PolarAxes = 'PolarAxesRepresentation'
    # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    extractSurface1Display.ScaleTransferFunction.Points = [
        1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    extractSurface1Display.OpacityTransferFunction.Points = [
        1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Hide data in view
    paraview.simple.Hide(shrink1, renderView1)
    # Update the view to ensure updated data information
    renderView1.Update()
    # create a new 'Tube'
    tube1 = paraview.simple.Tube(Input=extractSurface1)
    tube1.Scalars = [
        'POINTS', 'network |' + network.name +  '| labels | pore.all']
    tube1.Vectors = [None, '1']
    tube1.Radius = 0.04
    # Set active source
    paraview.simple.SetActiveSource(extractSurface1)
    # Destroy tube1
    paraview.simple.Delete(tube1)
    del tube1
    # Set active source
    paraview.simple.SetActiveSource(shrink1)
    # Set active source
    paraview.simple.SetActiveSource(cellDatatoPointData1)
    # Set active source
    paraview.simple.SetActiveSource(extractSurface1)
    # Create a new 'Tube'
    tube1 =paraview.simple.Tube(Input=extractSurface1)
    tube1.Scalars = [
        'POINTS', 'network |' + network.name + ' | labels | pore.all']
    tube1.Vectors = [None, '1']
    tube1.Radius = 0.04
    # Properties modified on tube1
    tube1.Vectors = ['POINTS', '1']
    # Show data in view
    tube1Display = paraview.simple.Show(tube1, renderView1,
                                        'GeometryRepresentation')
    # Trace defaults for the display properties.
    tube1Display.Representation = 'Surface'
    tube1Display.ColorArrayName = [None, '']
    tube1Display.OSPRayScaleArray = 'TubeNormals'
    tube1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tube1Display.SelectOrientationVectors = 'None'
    tube1Display.ScaleFactor = (maxshape)/10
    tube1Display.SelectScaleArray = 'None'
    tube1Display.GlyphType = 'Arrow'
    tube1Display.GlyphTableIndexArray = 'None'
    tube1Display.GaussianRadius = (maxshape)/200
    tube1Display.SetScaleArray = ['POINTS', 'TubeNormals']
    tube1Display.ScaleTransferFunction = 'PiecewiseFunction'
    tube1Display.OpacityArray = ['POINTS', 'TubeNormals']
    tube1Display.OpacityTransferFunction = 'PiecewiseFunction'
    tube1Display.DataAxesGrid = 'GridAxesRepresentation'
    tube1Display.PolarAxes = 'PolarAxesRepresentation'
    # Init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
    tube1Display.ScaleTransferFunction.Points = [-1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
    tube1Display.OpacityTransferFunction.Points = [-1, 0, 0.5, 0, 1, 1, 0.5, 0]
    # Hide data in view
    paraview.simple.Hide(extractSurface1, renderView1)
    # Update the view to ensure updated data information
    renderView1.Update()
    # Saving camera placements for all active views
    # Current camera placement for renderView1
    renderView1.CameraPosition = [
        (xshape + 1) / 2, (yshape + 1) / 2, 4.3 * np.sqrt(np.sum(shape / 2)**2)
    ]
    renderView1.CameraFocalPoint = [(xi+1) / 2 for xi in shape]
    renderView1.CameraParallelScale = np.sqrt(np.sum(shape / 2)**2)
    paraview.simple.SaveState(f"{file}.pvsm")


def open_paraview(filename):
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
