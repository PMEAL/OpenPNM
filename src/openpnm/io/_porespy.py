import pickle as pk
from openpnm.network import Network
from openpnm.io import _parse_filename


def network_from_porespy(filename):
    r"""
    Load a network extracted using the PoreSpy package

    Parameters
    ----------
    filename : str or dict
        Can either be a filename pointing to a pickled dictionary, or a
        handle to a dictionary in memory. The second option lets users
        avoid the step of saving the dictionary to a file.

    Returns
    -------
    network : dict
        An OpenPNM network dictionary

    """
    # Parse the filename
    if isinstance(filename, dict):
        net = dict(filename)
    else:
        filename = _parse_filename(filename=filename)
        with open(filename, mode='rb') as f:
            net = pk.load(f)

    # Pop param.* entries so they're routed via __setitem__ into _params
    # rather than landing as raw dict items via update().
    params = {k: net.pop(k) for k in list(net) if k.startswith('param.')}
    network = Network()
    network.update(net)
    for key, value in params.items():
        network[key] = value

    return network
