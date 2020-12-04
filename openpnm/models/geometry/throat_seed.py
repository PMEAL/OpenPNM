from openpnm.models import misc as _misc


def random(target, seed=None, num_range=[0, 1]):
    return _misc.random(target, element='throat', seed=seed,
                        num_range=num_range)


random.__doc__ = _misc.random.__doc__


def from_neighbor_pores(target, prop='pore.seed', mode='min'):
    return _misc.from_neighbor_pores(target=target, prop=prop,
                                     mode=mode)


random.__doc__ = _misc.from_neighbor_pores.__doc__
