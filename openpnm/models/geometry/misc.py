r"""
Helper methods for openpnm.geometry module.

"""


def _get_conduit_diameters(target, pore_diameter, throat_diameter):
    r"""
    Private helper methods for fetching conduit diameters.
    """
    network = target.project.network
    throats = target.throats(target=network)
    cn = network["throat.conns"][throats]

    D1 = network[pore_diameter][cn[:, 0]]
    Dt = network[throat_diameter][throats]
    D2 = network[pore_diameter][cn[:, 1]]

    return D1, Dt, D2
