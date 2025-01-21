import numpy as np
from uuid import uuid4
from openpnm.pnmlib.core import get_data, set_data, count


__all__ = [
    'add_network',
    'create_phase',
    'get_group',
]


def get_group(project, name):
    if name in ['/', '', '*']:
        return project
    temp = get_data(project, name+"/*")
    group = {}
    for k, v in temp.items():
        if not hasattr(v, 'keys'):
            group[k.rsplit('/', 1)[1]] = v
    return group


def add_network(project, network):
    coords = network.pop('pore.coords')
    set_data(project, 'pore.x', coords[:, 0])
    set_data(project, 'pore.y', coords[:, 1])
    set_data(project, 'pore.z', coords[:, 2])
    conns = network.pop('throat.conns')
    set_data(project, 'throat.t', conns[:, 0])
    set_data(project, 'throat.h', conns[:, 1])
    for k, v in network.items():
        if v.dtype == bool:
            set_data(project, k, v)
        else:
            set_data(project, 'network' + '/' + k, v)


def create_phase(target, name=None):
    if name is None:
        name = generate_name(target, prefix='phase')
    create_group(target=target, name=name)


def generate_name(target, prefix):
    i = 1
    while True:
        if prefix + '_' + str(i).zfill(2) not in target.keys():
            name = prefix + '_' + str(i).zfill(2)
            break
        i += 1
    return name


def create_group(target, name):
    for element in ['pore', 'throat']:
        try:
            n = count(target, element)
        except Exception:
            n = count(target['network'], element)
        set_data(target, name + '/' + element + '.all', np.ones((n, ), dtype=bool))
    return target
