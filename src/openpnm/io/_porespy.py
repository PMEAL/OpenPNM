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
        net = filename
    else:
        filename = _parse_filename(filename=filename)
        with open(filename, mode='rb') as f:
            net = pk.load(f)

    network = Network()
    network.update(net)

    return network
