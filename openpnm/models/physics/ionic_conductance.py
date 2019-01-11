from .misc import generic_conductance


def ordinary(target,
             pore_area='pore.area',
             throat_area='throat.area',
             conduit_lengths='throat.conduit_lengths'):
    r"""
    """
    return generic_conductance(target=target, transport_type='ionic',
                               pore_area=pore_area,
                               throat_area=throat_area,
                               pore_diffusivity='',
                               throat_diffusivity='',
                               conduit_lengths=conduit_lengths,
                               conduit_shape_factors='')
