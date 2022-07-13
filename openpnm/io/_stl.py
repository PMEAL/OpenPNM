import logging
import numpy as _np
from openpnm.io import _parse_filename
logger = logging.getLogger(__name__)


def network_to_stl(network, filename=None, maxsize='auto',
                   fileformat='STL Format', logger_level=0):
    r"""
    Saves (transient/steady-state) data from the given objects into the
    specified file.

    Parameters
    ----------
    network : Network.
        The network containing the desired data.
    phases : list[Phase] (place holder, default is None).
        List of phases containing the desired data.
    filename : str (optional).
        The name of the file containing the data to export.
    maxsize : a float or a string "auto" (optional).
        The maximum size of the mesh elements allowed. "auto" corresponds
        to an automatic determination based on pores and throats sizes. Any
        float value will be used as a maximum size. Small values result in
        finner meshes, but slower mesh calculations.
    fileformat : str (optional).
        Default is "STL Format" which corresponds to STL format. Other
        formats such as Gmsh and Fluent are supported (see ngsolve.org).
    logger_level : integer between 0 and 7 (optional).
        Default is 0. The logger level set in netgen package.

    Notes
    -----
    The STL Format is a Standard Triangle (or Tessellation) Language supported
    by many CAD packages and used for 3D printing. Export to this format may be
    slow since a 3D closed surface mesh is built.
    Requires installation of "netgen".
    """
    try:
        import netgen.csg as csg
    except ModuleNotFoundError:
        raise Exception('Module "netgen" not found.')
    try:
        from netgen.meshing import SetMessageImportance as log
        log(logger_level)
    except ModuleNotFoundError:
        logger.warning('Module "netgen.meshing" not found.'
                       + ' The "logger_level" will be ignored.')

    # Temporarily add endpoints to network so STL class works
    network["throat.endpoints.head"] = network.coords[network.conns[:, 0]]
    network["throat.endpoints.tail"] = network.coords[network.conns[:, 1]]

    filename = network.name if filename is None else filename
    path = _parse_filename(filename=filename, ext='stl')
    # Path is a pathlib object, so slice it up as needed
    fname_stl = path.name

    # Correct connections where 'pore.diameter' = 'throat.diameter'
    dt = network['throat.diameter'].copy()
    dp = network['pore.diameter'][network['throat.conns']]
    dt[dp[:, 0] == dt] *= 0.99
    dt[dp[:, 1] == dt] *= 0.99

    scale = max(network['pore.diameter'].max(), dt.max(),
                network['throat.length'].max())
    if maxsize == 'auto':
        maxsize = min(network['pore.diameter'].min(), dt.min(),
                      network['throat.length'].min())
    geo = csg.CSGeometry()

    # Define pores
    geometry = csg.Sphere(csg.Pnt(network['pore.coords'][0, 0]/scale,
                                  network['pore.coords'][0, 1]/scale,
                                  network['pore.coords'][0, 2]/scale),
                          network['pore.diameter'][0]/scale/2)
    for p in range(1, network.Np):
        pore = csg.Sphere(csg.Pnt(network['pore.coords'][p, 0]/scale,
                                  network['pore.coords'][p, 1]/scale,
                                  network['pore.coords'][p, 2]/scale),
                          network['pore.diameter'][p]/scale/2)
        geometry += pore

    # Define throats
    for t in range(network.Nt):
        A = network['throat.endpoints.tail'][t, :]/scale
        B = network['throat.endpoints.head'][t, :]/scale
        V = (B-A)/_np.linalg.norm(B-A)
        plane1 = csg.Plane(csg.Pnt(A[0], A[1], A[2]),
                           csg.Vec(-V[0], -V[1], -V[2]))
        plane2 = csg.Plane(csg.Pnt(B[0], B[1], B[2]),
                           csg.Vec(V[0], V[1], V[2]))
        cylinder = csg.Cylinder(csg.Pnt(A[0], A[1], A[2]),
                                csg.Pnt(B[0], B[1], B[2]),
                                dt[t]/scale/2)
        throat = cylinder * plane1 * plane2
        geometry += throat

    # Add pore and throats to geometry, build mesh, rescale, and export
    geo.Add(geometry)
    mesh = geo.GenerateMesh(maxh=maxsize/scale)
    mesh.Scale(scale)
    mesh.Export(filename=fname_stl, format=fileformat)

    # Remove endpoints label from network
    del network["throat.endpoints"]
