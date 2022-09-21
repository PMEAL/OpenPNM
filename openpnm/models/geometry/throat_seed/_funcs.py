from openpnm.models import misc as _misc


__all__ = ["random",
           "from_neighbor_pores"]


def random(network, seed=None, num_range=[0, 1]):
    return _misc.random(network, element='throat', seed=seed,
                        num_range=num_range)


random.__doc__ = _misc.random.__doc__


def from_neighbor_pores(network, prop='pore.seed', mode='min'):
    return _misc.from_neighbor_pores(network, prop=prop,
                                     mode=mode)


random.__doc__ = _misc.from_neighbor_pores.__doc__
