import numpy as np
try:
    import vtk
except ImportError:
    pass

def preview(pn, values=[]):

    coords = pn.get_pore_data(prop='coords')
    heads, tails = pn.get_throat_data(prop='conns').T

    points = vtk.vtkPoints()
    for x,y,z in coords:
        points.InsertNextPoint(x, y, z)

    polys = vtk.vtkCellArray()
    for hi, ti in zip(heads, tails):
        vil = vtk.vtkIdList()
        vil.InsertNextId(hi)
        vil.InsertNextId(ti)
        polys.InsertNextCell(vil)

    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    if any(values):
        values = np.array(values)
        values -= values.min()
        values = np.true_divide(values, values.max())
    for v in values:
        r = 255*(v)
        g = 0
        b = 255*(1-v)
        colors.InsertNextTuple3(r,g,b)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(polys)
    if colors.GetNumberOfTuples() == len(coords):
        polydata.GetPointData().SetScalars(colors)
    elif colors.GetNumberOfTuples() != 0:
        raise Exception("Mismatch: {} points, {} scalars".format(
                        len(coords), colors.GetNumberOfTuples()))

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInput(polydata)

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    ren = vtk.vtkRenderer()
    ren.AddActor(actor)

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    iren.Initialize()
    renWin.Render()
    iren.Start()