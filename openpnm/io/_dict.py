import logging
from openpnm.utils import NestedDict


logger = logging.getLogger(__name__)


def project_to_dict(project, categorize_by=[], flatten=False, element=None,
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

    if flatten:
        d = {}
    else:
        d = NestedDict(delimiter=delim)

    def build_path(obj, key):
        propname = delim + key
        prefix = ''
        datatype = ''
        arr = obj[key]
        if 'object' in categorize_by:
            if hasattr(obj, 'coords'):
                prefix = 'network' + delim
            else:
                prefix = 'phase' + delim
        if 'element' in categorize_by:
            propname = delim + key.replace('.', delim)
        if 'data' in categorize_by:
            if arr.dtype == bool:
                datatype = delim + 'labels'
            else:
                datatype = delim + 'properties'
        path = prefix + obj.name + datatype + propname
        return path

    for net in network:
        for key in net.props(element=element) + net.labels(element=element):
            path = build_path(obj=net, key=key)
            d[path] = net[key]

    for phase in phases:
        for key in phase.props(element=element) + phase.labels(element=element):
            path = build_path(obj=phase, key=key)
            d[path] = phase[key]

    return d
