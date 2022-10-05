import numpy as np


def cubic(shape, spacing=1, connectivity=6, node_prefix='node', edge_prefix='edge'):
    r"""
    Generate a simple cubic lattice

    Parameters
    ----------
    shape : array_like
        The number of unit cells in each direction.  A unit cell has 1 vertex
        at its center.
    spacing : array_like or float
        The size of a unit cell in each direction. If an scalar is given it is
        applied in all 3 directions.

    Returns
    -------
    network : dict
        A dictionary containing ``coords`` and ``conns`` of a cubic network with the
        specified spacing and connectivity.

    """
    # Take care of 1D/2D networks
    shape = np.array(shape, ndmin=1)
    shape = np.concatenate((shape, [1] * (3 - shape.size))).astype(int)
    arr = np.atleast_3d(np.empty(shape))
    spacing = np.float64(spacing)
    if spacing.size == 2:
        spacing = np.concatenate((spacing, [1]))
    spacing = np.ones(3, dtype=float) * np.array(spacing, ndmin=1)

    z = np.tile(np.arange(shape[2]), shape[0] * shape[1])
    y = np.tile(np.repeat(np.arange(shape[1]), shape[2]), shape[0])
    x = np.repeat(np.arange(shape[0]), shape[1] * shape[2])
    points = (np.vstack([x, y, z]).T).astype(float) + 0.5

    idx = np.arange(arr.size).reshape(arr.shape)

    face_joints = [(idx[:, :, :-1], idx[:, :, 1:]),
                   (idx[:, :-1], idx[:, 1:]),
                   (idx[:-1], idx[1:])]

    corner_joints = [(idx[:-1, :-1, :-1], idx[1:, 1:, 1:]),
                     (idx[:-1, :-1, 1:], idx[1:, 1:, :-1]),
                     (idx[:-1, 1:, :-1], idx[1:, :-1, 1:]),
                     (idx[1:, :-1, :-1], idx[:-1, 1:, 1:])]

    edge_joints = [(idx[:, :-1, :-1], idx[:, 1:, 1:]),
                   (idx[:, :-1, 1:], idx[:, 1:, :-1]),
                   (idx[:-1, :, :-1], idx[1:, :, 1:]),
                   (idx[1:, :, :-1], idx[:-1, :, 1:]),
                   (idx[1:, 1:, :], idx[:-1, :-1, :]),
                   (idx[1:, :-1, :], idx[:-1, 1:, :])]

    if connectivity == 6:
        joints = face_joints
    elif connectivity == 6 + 8:
        joints = face_joints + corner_joints
    elif connectivity == 6 + 12:
        joints = face_joints + edge_joints
    elif connectivity == 12 + 8:
        joints = edge_joints + corner_joints
    elif connectivity == 6 + 8 + 12:
        joints = face_joints + corner_joints + edge_joints
    else:
        raise Exception("Invalid connectivity. Must be 6, 14, 18, 20 or 26.")

    tails, heads = np.array([], dtype=int), np.array([], dtype=int)
    for T, H in joints:
        tails = np.concatenate((tails, T.flatten()))
        heads = np.concatenate((heads, H.flatten()))
    pairs = np.vstack([tails, heads]).T
    # NOTE: pairs is already sorted for connectivity = 6
    if connectivity != 6:
        pairs = np.sort(pairs, axis=1)

    d = {}
    d[f"{node_prefix}.coords"] = points * spacing
    d[f"{edge_prefix}.conns"] = pairs

    return d
