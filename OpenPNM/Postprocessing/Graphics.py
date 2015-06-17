from tempfile import NamedTemporaryFile
from subprocess import call
from itertools import count
import numpy as np
from matplotlib import cm

try:
    import vtk
except ImportError:
    vtk = type('', (), {'vtkActor': object})


class Actor(vtk.vtkActor):

    def __init__(self):
        raise NotImplementedError()

    def pointArray(self, _points):
        points = vtk.vtkPoints()
        for x, y, z in _points:
            points.InsertNextPoint(x, y, z)
        return points

    def lineArray(self, _lines):
        lines = vtk.vtkCellArray()
        for id_set in _lines:
            l = vtk.vtkIdList()
            for i in id_set:
                l.InsertNextId(i)
            lines.InsertNextCell(l)
        return lines

    def faceArray(self, faces):
        lines = vtk.vtkCellArray()
        for face in faces:
            l = vtk.vtkIdList()
            for i in face:
                l.InsertNextId(i)
            lines.InsertNextCell(l)
        return lines

    def floatArray(self, array):
        floats = vtk.vtkFloatArray()
        for n in array:
            floats.InsertNextValue(n)
        return floats

    def colorArray(self, array, cmap=None):
        colors = vtk.vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        if cmap is None:
            cmap = 'coolwarm'
        colormap = cm.get_cmap(cmap)
        mapped = colormap(array)
        for r, g, b, a in 255*mapped:
            colors.InsertNextTuple3(r, g, b)
        return colors

    def update(self, t=0):
        pass


class Wires(Actor):

    def __init__(self, vertex_coords, edge_pairs, vertex_weights=None,
                 alpha=1, cmap=None):
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(self.pointArray(vertex_coords))
        self.polydata.SetLines(self.lineArray(edge_pairs))
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.polydata)
        self.SetMapper(self.mapper)

        if vertex_weights is None:
            weights = np.atleast_2d([128 for _ in vertex_coords])
        else:
            weights = np.atleast_2d(vertex_weights)
            weights = np.subtract(weights, weights.min())
            weights = np.true_divide(weights, weights.max())
        self.weights = weights
        self.cmap = cmap
        self.update()

        self.GetProperty().SetOpacity(alpha)

    def update(self, t=0):
        i = t % len(self.weights)
        self.polydata.GetPointData().SetScalars(self.colorArray(self.weights[i],
                                                                self.cmap))


class Surface(Actor):

    def __init__(self, vertex_coords, edge_pairs, vertex_weights=None, alpha=1,
                 cmap=None, offset=0):
        self.offset = offset
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(self.pointArray(vertex_coords))
        self.polydata.SetLines(self.lineArray(edge_pairs))
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInput(self.polydata)
        self.SetMapper(self.mapper)

        if vertex_weights is None:
            weights = np.atleast_2d([128 for _ in vertex_coords])
        else:
            weights = np.atleast_2d(vertex_weights)
        self.weights = weights
        self.cmap = cmap
        self.update()

        self.GetProperty().SetOpacity(alpha)

    def update(self, t=0):
        i = t % len(self.weights)
        values = self.weights[i]
        # update z-positions
        points = self.polydata.GetPoints()
        for i, v in enumerate(values):
            x, y, _ = points.GetPoint(i)
            points.SetPoint(i, x, y, v + self.offset)
        # update colors
        normalized = np.true_divide(values, np.abs(self.weights).max())
        normalized = np.subtract(normalized, self.weights.min())
        colors = self.colorArray(normalized, self.cmap)
        self.polydata.GetPointData().SetScalars(colors)


class Tubes(Actor):

    def __init__(self, centers, vectors, radii, alpha=1, cmap=None):
        tails = centers - vectors/2.
        heads = centers + vectors/2.
        points = np.vstack(zip(tails, heads))
        pairs = np.arange(len(centers)*2).reshape(-1, 2)
        radii = radii.repeat(2)

        assert (points.size/3. == pairs.size)
        assert (pairs.size == radii.size)

        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(self.pointArray(points))
        self.polydata.SetLines(self.lineArray(pairs))
        self.polydata.GetPointData().SetScalars(self.floatArray(radii))

        self.tubeFilter = vtk.vtkTubeFilter()
        self.tubeFilter.SetInput(self.polydata)
        self.tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
        self.tubeFilter.SetNumberOfSides(10)
        self.tubeFilter.CappingOn()

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.tubeFilter.GetOutputPort())
        self.mapper.ScalarVisibilityOff()
        self.SetMapper(self.mapper)

        self.GetProperty().SetOpacity(alpha)


class Spheres(Actor):

    def __init__(self, centers, radii, alpha=1, color=(1, 1, 1)):
        self.polydata = vtk.vtkPolyData()
        self.polydata.SetPoints(self.pointArray(centers))
        self.radii = np.atleast_2d(radii)
        self.update()

        self.sphere_source = vtk.vtkSphereSource()
        self.glypher = vtk.vtkProgrammableGlyphFilter()
        self.glypher.SetInput(self.polydata)
        self.glypher.SetSource(self.sphere_source.GetOutput())
        self.glypher.SetGlyphMethod(self.glyph_method)

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.glypher.GetOutputPort())
        self.SetMapper(self.mapper)

        self.GetProperty().SetOpacity(alpha)
        r, g, b = color
        self.mapper.SetScalarVisibility(False)
        self.GetProperty().SetColor(r, g, b)

    def glyph_method(self):
        pid = self.glypher.GetPointId()
        self.sphere_source.SetCenter(self.glypher.GetPoint())
        radius = self.glypher.GetPointData().GetScalars().GetValue(pid)
        self.sphere_source.SetRadius(radius)

    def update(self, t=0):
        i = t % len(self.radii)
        self.polydata.GetPointData().SetScalars(self.floatArray(self.radii[i]))


class Scene(object):
    ticks = count(0)

    def __init__(self, parent=None, fix_camera=True):
        """
        fix_camera : more sensible default
        """
        if parent is not None:
            self.renWin = parent.GetRenderWindow()
            self.iren = self.renWin.GetInteractor()
        else:
            self.renWin = vtk.vtkRenderWindow()
            self.iren = vtk.vtkRenderWindowInteractor()
            self.iren.SetRenderWindow(self.renWin)

        self.ren = vtk.vtkRenderer()
        self.renWin.AddRenderer(self.ren)

        if fix_camera:
            camera = vtk.vtkInteractorStyleTrackballCamera()
            self.iren.SetInteractorStyle(camera)

    def update_all(self, object=None, event=None, t=None):
        if t is None:
            t = next(self.ticks)
        for aid in range(self.ren.VisibleActorCount()):
            actor = self.ren.GetActors().GetItemAsObject(aid)
            actor.update(t)
        self.renWin.Render()

    def save(self, frames, outfile='animated.gif'):
        """
        takes a snapshot of the frames at given t, and returns the paths
        """
        windowToImage = vtk.vtkWindowToImageFilter()
        windowToImage.SetInput(self.renWin)
        writer = vtk.vtkPNGWriter()
        writer.SetInput(windowToImage.GetOutput())

        slide_paths = []
        for t in frames:
            f = NamedTemporaryFile(suffix='.png', delete=False)
            self.update_all(t=t)
            windowToImage.Modified()
            writer.SetFileName(f.name)
            writer.Write()
            slide_paths.append(f.name)

        call(["convert"] + slide_paths + [outfile])
        call(["rm"] + slide_paths)

    def play(self, timeout=1):
        self.iren.Initialize()
        if timeout is not None:
            self.iren.AddObserver('TimerEvent', self.update_all)
            self.timer = self.iren.CreateRepeatingTimer(timeout)
        self.update_all()
        self.iren.Start()

    def add_actors(self, list_of_actors):
        for actor in list_of_actors:
            self.ren.AddActor(actor)
