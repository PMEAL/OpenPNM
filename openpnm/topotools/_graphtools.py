from openpnm._skgraph import tools
from openpnm._skgraph import queries
from openpnm._skgraph import operations


__all__ = [
    'find_neighbor_sites',
    'find_neighbor_bonds',
    'find_connected_sites',
    'find_connecting_bonds',
    'istriu',
    'istril',
    'istriangular',
    'issymmetric',
    'tri_to_am',
    'vor_to_am',
    'conns_to_am',
    'drop_sites',
]


def find_neighbor_sites(sites, **kwargs):
    return queries.find_neighbor_nodes(inds=sites, **kwargs)


find_neighbor_sites.__doc__ = queries.find_neighbor_nodes.__doc__


def find_neighbor_bonds(sites, **kwargs):
    return queries.find_neighbor_edges(inds=sites, **kwargs)


find_neighbor_bonds.__doc__ = queries.find_neighbor_edges.__doc__


def find_connected_sites(bonds, **kwargs):
    return queries.find_connected_nodes(inds=bonds, **kwargs)


find_connected_sites.__doc__ = queries.find_connected_nodes.__doc__


def find_connecting_bonds(sites, **kwargs):
    return queries.find_connecting_edges(inds=sites, **kwargs)


find_connecting_bonds.__doc__ = queries.find_connecting_edges.__doc__


def istriu(am):
    return tools.istriu(am)


istriu.__doc__ = tools.istriu.__doc__


def istril(am):
    return tools.istril(am)


istril.__doc__ = tools.istril.__doc__


def istriangular(am):
    return tools.istriangular(am)


istriangular.__doc__ = tools.istriangular.__doc__


def issymmetric(am):
    return tools.issymmetric(am)


issymmetric.__doc__ = tools.issymmetric.__doc__


def tri_to_am(tri):
    return tools.tri_to_am(tri=tri)


tri_to_am.__doc__ = tools.tri_to_am.__doc__


def vor_to_am(vor):
    return tools.vor_to_am(vor=vor)


vor_to_am.__doc__ = tools.vor_to_am.__doc__


def conns_to_am(*args, **kwargs):
    return tools.conns_to_am(*args, **kwargs)


conns_to_am.__doc__ = tools.conns_to_am.__doc__


def drop_sites(am, sites):
    return operations.drop_nodes_from_am(inds=sites, am=am)


drop_sites.__doc__ = operations.drop_nodes_from_am.__doc__
