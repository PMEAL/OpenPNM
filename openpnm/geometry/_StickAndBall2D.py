import openpnm.models as mods
from openpnm.geometry import GenericGeometry
from numpy.linalg import norm as _norm
from openpnm.models.geometry.throat_length import ctc as _ctc
import numpy as _np



class _StickAndBall2D(GenericGeometry):
    r"""
    2D Stick and Ball subclass of GenericGeometry.  This subclass is meant as a
    basic default geometry to get started quickly.

    Pore diameters are randomly assigned between 0 and the largest sphere that
    does not overlap with it's nearest neighbor.

    Throat diameters are half the diameter of the smaller of it's two
    neighboring pores.

    Parameters
    ----------
    network : OpenPNM Network object
        The network with which this Geometry should be associated

    project : OpenPNM Project object, optional
        Can be supplied instead of a ``network``

    pores : array_like
        The pores in the domain where this Geometry applies

    throats : array_like
        The throats in the domain where this Geometry applies

    name : string
        The name of the object, which is also used as the label where this
        geometry is defined.

    Examples
    --------
    Geometry objects (along with Physics objects) can be applied to a subset
    of pores and/or throats.  This allows for different geometrical property
    models to be applied in different regions.  This is illustrated in the
    following code:

    >>> import numpy as np
    >>> import scipy as sp
    >>> import openpnm as op
    >>> import matplotlib.pyplot as plt
    >>> pn = op.network.CubicDual(shape=[5, 5, 5])
    >>> Ps = pn.pores('primary')
    >>> Ts = pn.throats('primary')
    >>> geo1 = op.geometry._StickAndBall2D(network=pn, pores=Ps, throats=Ts)
    >>> Ps = pn.pores('secondary')
    >>> Ts = pn.throats(['secondary', 'interconnect'])
    >>> geo2 = op.geometry._StickAndBall2D(network=pn, pores=Ps, throats=Ts)

    Now override the 'pore.diameter' values on the ``geo2`` object:

    >>> geo2.remove_model('pore.diameter')  # Remove model and data
    >>> geo2['pore.diameter'] = np.random.rand(geo2.Np)*0.05

    Look at the 'pore.diameter' distributions on each object:

    >>> fig = plt.hist(geo1['pore.diameter'], bins=20, alpha=0.5)
    >>> fig = plt.hist(geo2['pore.diameter'], bins=20, alpha=0.5)

    The resulting figure shows that these two Geometry object each have a
    different pore size distribution, with ``geo2`` being much smaller:

    .. image:: /../docs/_static/images/stick_and_ball_histogram.png
        :align: center

    Notes
    -----
    The table below gives a summary of the all the pore-scale models that are
    included on this class.

    All of these parameters can be adjusted manually by editing the entries in
    the **ModelsDict** stored in the ``models`` attribute of the object.

    For a full listing of models and their parameters use ``print(obj.models)``
    where ``obj`` is the handle to the object.

    In addition to these models, this class also has a number of constant
    values assigned to it which can be found by running
    ``props(mode='constants')``.

    +----+----------------------+------------------+--------------------------+
    | #  | Property Name        | Parameter        | Value                    |
    +====+======================+==================+==========================+
    | 1  | pore.seed            | model:           | random                   |
    +----+----------------------+------------------+--------------------------+
    |    |                      | element          | pore                     |
    +----+----------------------+------------------+--------------------------+
    |    |                      | num_range        | [0, 0.1]                 |
    +----+----------------------+------------------+--------------------------+
    |    |                      | seed             | None                     |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 2  | pore.max_size        | model:           | largest_sphere           |
    +----+----------------------+------------------+--------------------------+
    |    |                      | iters            | 10                       |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    |    |                      | fixed_diameter   | pore.fixed_diameter      |
    +----+----------------------+------------------+--------------------------+
    | 3  | pore.diameter        | model:           | product                  |
    +----+----------------------+------------------+--------------------------+
    |    |                      | prop1            | pore.max_size            |
    +----+----------------------+------------------+--------------------------+
    |    |                      | prop2            | pore.seed                |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 4  | pore.area            | model:           | sphere                   |
    +----+----------------------+------------------+--------------------------+
    |    |                      | pore_diameter    | pore.diameter            |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 5  | pore.volume          | model:           | sphere                   |
    +----+----------------------+------------------+--------------------------+
    |    |                      | pore_diameter    | pore.diameter            |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 6  | throat.max_size      | model:           | from_neighbor_pores      |
    +----+----------------------+------------------+--------------------------+
    |    |                      | mode             | min                      |
    +----+----------------------+------------------+--------------------------+
    |    |                      | pore_prop        | pore.diameter            |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 7  | throat.diameter      | model:           | scaled                   |
    +----+----------------------+------------------+--------------------------+
    |    |                      | factor           | 0.5                      |
    +----+----------------------+------------------+--------------------------+
    |    |                      | prop             | throat.max_size          |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 8  | throat.length        | model:           | piecewise                |
    +----+----------------------+------------------+--------------------------+
    |    |                      | pore_diameter    | pore.diameter            |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 9  | throat.surface_area  | model:           | cylinder                 |
    +----+----------------------+------------------+--------------------------+
    |    |                      | throat_diameter  | throat.diameter          |
    +----+----------------------+------------------+--------------------------+
    |    |                      | throat_length    | throat.length            |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 10 | throat.volume        | model:           | cylinder                 |
    +----+----------------------+------------------+--------------------------+
    |    |                      | throat_diameter  | throat.diameter          |
    +----+----------------------+------------------+--------------------------+
    |    |                      | throat_length    | throat.length            |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+
    | 11 | throat.area          | model:           | cylinder                 |
    +----+----------------------+------------------+--------------------------+
    |    |                      | throat_diameter  | throat.diameter          |
    +----+----------------------+------------------+--------------------------+
    |    |                      | regen_mode       | normal                   |
    +----+----------------------+------------------+--------------------------+

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self.add_model(propname='pore.seed',
                       model=mods.misc.random,
                       element='pore',
                       num_range=[0.2, 0.7],
                       seed=None)

        self.add_model(propname='pore.max_size',
                       model=mods.geometry.pore_size.largest_sphere,
                       iters=10)

        self.add_model(propname='pore.diameter',
                       model=mods.misc.product,
                       props=['pore.max_size', 'pore.seed'])

        self.add_model(propname='pore.area',
                       model=mods.geometry.pore_cross_sectional_area.circle,
                       pore_diameter='pore.diameter')

        self.add_model(propname='pore.volume',
                       model=mods.geometry.pore_volume.circle,
                       pore_diameter='pore.diameter')

        self.add_model(propname='throat.max_size',
                       model=mods.misc.from_neighbor_pores,
                       mode='min',
                       prop='pore.diameter')

        self.add_model(propname='throat.diameter',
                       model=mods.misc.scaled,
                       factor=0.5,
                       prop='throat.max_size')

        self.add_model(propname='throat.endpoints',
                       model=circular_pores,
                       pore_diameter='pore.diameter',
                       throat_diameter='throat.diameter')

        self.add_model(propname='throat.length',
                       model=piecewise,
                       throat_endpoints='throat.endpoints')

        self.add_model(propname='throat.surface_area',
                       model=mods.geometry.throat_surface_area.rectangle,
                       throat_length='throat.length')

        self.add_model(propname='throat.volume',
                       model=mods.geometry.throat_volume.rectangle,
                       throat_diameter='throat.diameter',
                       throat_length='throat.length')

        self.add_model(propname='throat.area',
                       model=mods.geometry.throat_cross_sectional_area.rectangle,
                       throat_diameter='throat.diameter')

        self.add_model(propname='throat.conduit_lengths',
                       model=conduit_lengths,
                       throat_endpoints='throat.endpoints',
                       throat_length='throat.length')

def circular_pores(target, pore_diameter='pore.diameter',
                    throat_diameter='throat.diameter',
                    throat_centroid='throat.centroid'):
    r"""
    Calculate the coordinates of throat endpoints, assuming spherical pores.
    This model accounts for the overlapping lens between pores and throats.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_diameter : string
        Dictionary key of the pore diameter values.

    throat_diameter : string
        Dictionary key of the throat diameter values.

    throat_centroid : string, optional
        Dictionary key of the throat centroid values. See the notes.

    Returns
    -------
    EP : dictionary
        Coordinates of throat endpoints stored in Dict form. Can be accessed
        via the dict keys 'head' and 'tail'.

    Notes
    -----
    (1) This model should not be applied to true 2D networks. Use
    ``circular_pores`` model instead.

    (2) By default, this model assumes that throat centroid and pore
    coordinates are colinear. If that's not the case, such as in extracted
    networks, ``throat_centroid`` could be passed as an optional argument, and
    the model takes care of the rest.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    xyz = network['pore.coords']
    cn = network['throat.conns'][throats]
    L = _ctc(target=target) + 1e-15
    Dt = network[throat_diameter][throats]
    D1 = network[pore_diameter][cn[:, 0]]
    D2 = network[pore_diameter][cn[:, 1]]
    L1 = _np.zeros_like(L)
    L2 = _np.zeros_like(L)
    # Handle the case where Dt > Dp
    mask = Dt > D1
    L1[mask] = 0.5 * D1[mask]
    L1[~mask] = _np.sqrt(D1[~mask]**2 - Dt[~mask]**2) / 2
    mask = Dt > D2
    L2[mask] = 0.5 * D2[mask]
    L2[~mask] = _np.sqrt(D2[~mask]**2 - Dt[~mask]**2) / 2
    # Handle non-colinear pores and throat centroids
    try:
        TC = network[throat_centroid][throats]
        LP1T = _np.linalg.norm(TC - xyz[cn[:, 0]], axis=1) + 1e-15
        LP2T = _np.linalg.norm(TC - xyz[cn[:, 1]], axis=1) + 1e-15
        unit_vec_P1T = (TC - xyz[cn[:, 0]]) / LP1T[:, None]
        unit_vec_P2T = (TC - xyz[cn[:, 1]]) / LP2T[:, None]
    except KeyError:
        unit_vec_P1T = (xyz[cn[:, 1]] - xyz[cn[:, 0]]) / L[:, None]
        unit_vec_P2T = -1 * unit_vec_P1T
    # Find throat endpoints
    EP1 = xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T
    EP2 = xyz[cn[:, 1]] + L2[:, None] * unit_vec_P2T
    # Handle throats w/ overlapping pores
    L1 = (4 * L**2 + D1**2 - D2**2) / (8 * L)
    L2 = (4 * L**2 + D2**2 - D1**2) / (8 * L)
    h = (2 * _np.sqrt(D1**2 / 4 - L1**2)).real
    overlap = L - 0.5 * (D1 + D2) < 0
    mask = overlap & (Dt < h)
    EP1[mask] = (xyz[cn[:, 0]] + L1[:, None] * unit_vec_P1T)[mask]
    EP2[mask] = (xyz[cn[:, 1]] + L2[:, None] * unit_vec_P2T)[mask]
    return {'head': EP1, 'tail': EP2}

def piecewise(
    target,
    throat_endpoints='throat.endpoints',
    throat_centroid='throat.centroid'
):
    r"""
    Calculates throat length from end points and optionally a centroid

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    throat_endpoints : str
        Dictionary key of the throat endpoint values.
    throat_centroid : str
        Dictionary key of the throat centroid values, optional.

    Returns
    -------
    Lt : ndarray
        Array containing throat lengths for the given geometry.

    Notes
    -----
    By default, the model assumes that the centroids of pores and the
    connecting throat in each conduit are colinear.

    If `throat_centroid` is passed, the model accounts for the extra
    length. This could be useful for Voronoi or extracted networks.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    # Get throat endpoints
    EP1 = network[throat_endpoints + '.head'][throats]
    EP2 = network[throat_endpoints + '.tail'][throats]
    # Calculate throat length
    Lt = _norm(EP1 - EP2, axis=1)
    # Handle the case where pores & throat centroids are not colinear
    try:
        Ct = network[throat_centroid][throats]
        Lt = _norm(Ct - EP1, axis=1) + _norm(Ct - EP2, axis=1)
    except KeyError:
        pass
    return Lt


def conduit_lengths(
    target,
    throat_endpoints='throat.endpoints',
    throat_length='throat.length',
    throat_centroid='throat.centroid'
):
    r"""
    Calculates conduit lengths. A conduit is defined as half pore + throat
    + half pore.

    Parameters
    ----------
    target : GenericGeometry
        Geometry object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    throat_endpoints : str
        Dictionary key of the throat endpoint values.

    throat_diameter : str
        Dictionary key of the throat length values.

    throat_length : string (optional)
        Dictionary key of the throat length values.  If not given then the
        direct distance bewteen the two throat end points is used.

    Returns
    -------
    Dictionary containing conduit lengths, which can be accessed via the dict
    keys 'pore1', 'pore2', and 'throat'.

    """
    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    cn = network['throat.conns'][throats]
    # Get pore coordinates
    C1 = network['pore.coords'][cn[:, 0]]
    C2 = network['pore.coords'][cn[:, 1]]
    # Get throat endpoints and length
    EP1 = network[throat_endpoints + '.head'][throats]
    EP2 = network[throat_endpoints + '.tail'][throats]
    try:
        # Look up throat length if given
        Lt = network[throat_length][throats]
    except KeyError:
        # Calculate throat length otherwise based on piecewise model
        Lt = piecewise(target, throat_endpoints, throat_centroid)
    # Calculate conduit lengths for pore 1 and pore 2
    L1 = _norm(C1 - EP1, axis=1)
    L2 = _norm(C2 - EP2, axis=1)
    Lt = _np.maximum(Lt, 1e-15)

    return {'pore1': L1, 'throat': Lt, 'pore2': L2}
