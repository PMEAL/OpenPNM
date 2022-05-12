__all__ = [
    'from_porespy',
    'from_marock',
    'from_networkx',
    'from_pergeos',
    'from_jsongraph',
    'from_csv',
]


def from_porespy(file):
    from openpnm.io import network_from_porespy
    net = network_from_porespy(file)
    return net


def from_marock(*args, **kwargs):
    from openpnm.io import network_from_marock
    net = network_from_marock(*args, **kwargs)
    return net


def from_networkx(*args, **kwargs):
    from openpnm.io import network_from_networkx
    net = network_from_networkx(*args, **kwargs)
    return net


def from_pergeos(*args, **kwargs):
    from openpnm.io import network_from_pergeos
    net = network_from_pergeos(*args, **kwargs)
    return net


def from_jsongraph(*args, **kwargs):
    from openpnm.io import network_from_jsongraph
    net = network_from_jsongraph(*args, **kwargs)
    return net


def from_csv(*args, **kwargs):
    from openpnm.io import network_from_csv
    net = network_from_csv(*args, **kwargs)
    return net
