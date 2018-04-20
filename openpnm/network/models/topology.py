import scipy as sp
import scipy.sparse.csgraph as csg


def coordination_number(target):
    Nn = target.num_neighbors(pores=target.Ps, flatten=False)
    return Nn


def reduce_coordination(target, z):
    # Find minimum spanning tree using random weights
    am = target.create_adjacency_matrix(sym=False, data=sp.rand(target.Nt))
    mst = csg.minimum_spanning_tree(am, overwrite=True)
    mst = mst.tocoo()
    # Label throats on spanning tree to avoid deleting them
    Ts = target.find_connecting_throat(mst.row, mst.col)
    Ts = sp.hstack(Ts)
    target['throat.mst'] = False
    target['throat.mst'][Ts] = True
    # Trim throats not on the spanning tree to acheive desired coordination
    Ts = sp.random.permutation(target.throats('mst', mode='complement'))
    Ts = Ts[:int(target.Nt - target.Np*(z/2))]
    mask = target.tomask(throats=Ts)
    return mask
