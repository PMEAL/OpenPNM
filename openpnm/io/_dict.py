import logging
from openpnm.utils import NestedDict
from openpnm.utils._misc import is_transient, nbr_to_str


logger = logging.getLogger(__name__)


def project_to_dict(project, categorize_by=['name'], flatten=False, element=None,
                    delim=' | '):
    r"""
    Returns a single dictionary object containing data from the given
    OpenPNM project, with the keys organized differently depending on
    optional arguments.

    Parameters
    ----------
    project : list
        An OpenPNM project object
    categorize_by : str or list[str]
        Indicates how the dictionaries should be organized.  The list can
        contain any, all or none of the following strings:

        **'object'** : If specified the dictionary keys will be stored
        under a general level corresponding to their type (e.g.
        'network/net_01/pore.all').

        **'name'** : If specified, then the data arrays are additionally
        categorized by their name.  This is enabled by default.

        **'data'** : If specified the data arrays are additionally
        categorized by ``label`` and ``property`` to separate *boolean*
        from *numeric* data.

        **'element'** : If specified the data arrays are
        additionally categorized by ``pore`` and ``throat``, meaning
        that the propnames are no longer prepended by a 'pore.' or
        'throat.'

    Returns
    -------
    A dictionary with the data stored in a hierarchical data structure, the
    actual format of which depends on the arguments to the function.

    """
    network = project.network
    phases = project.phases
    algs = project.algorithms

    if flatten:
        d = {}
    else:
        d = NestedDict(delimiter=delim)

    def build_path(obj, key):
        propname = key
        name = ''
        prefix = ''
        datatype = ''
        arr = obj[key]
        if 'object' in categorize_by:
            if hasattr(obj, 'coords'):
                prefix = 'network' + delim
            else:
                prefix = 'phase' + delim
        if 'element' in categorize_by:
            propname = key.replace('.', delim)
        if 'data' in categorize_by:
            if arr.dtype == bool:
                datatype = 'labels' + delim
            else:
                datatype = 'properties' + delim
        if 'name' in categorize_by:
            name = obj.name + delim
        path = prefix + name + datatype + propname
        return path

    for key in network.props(element=element) + network.labels(element=element):
        path = build_path(obj=network, key=key)
        d[path] = network[key]

    for phase in phases:
        for key in phase.props(element=element) + phase.labels(element=element):
            path = build_path(obj=phase, key=key)
            d[path] = phase[key]

    for alg in algs:
        try:  # 'quantity' is missing for multiphysics algorithm
            key = alg.settings['quantity']
        except AttributeError:
            break
        # only for transient algs
        transient = is_transient(alg)
        path = build_path(alg, key)
        if transient:
            times = alg.soln[key].t
            for i, t in enumerate(times):
                d[path + '#' + nbr_to_str(t)] = alg.soln[key][:, i]
        else:
            d[path] = alg.soln[key]

    return d
