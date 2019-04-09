r"""
.. autofunction:: openpnm.models.physics.diffusive_conductance.ordinary_diffusion
"""

from .misc import generic_conductance


def ad_dif_mig(target,
               pore_area='pore.area',
               throat_area='throat.area',
               pore_diffusivity='pore.diffusivity',
               throat_diffusivity='throat.diffusivity',
               conduit_lengths='throat.conduit_lengths',
               conduit_shape_factors='throat.poisson_shape_factors',
               pore_pressure='pore.pressure',
               pore_potential='pore.potential',
               throat_hydraulic_conductance='throat.hydraulic_conductance',
               throat_diffusive_conductance='throat.diffusive_conductance',
               throat_valence='throat.valence',
               throat_temperature='throat.temperature',
               electrolyte='',
               s_scheme='powerlaw'):
    r"""
    """
    return generic_conductance(
        target=target,
        transport_type='ad_dif_mig',
        pore_area=pore_area,
        throat_area=throat_area,
        pore_diffusivity=pore_diffusivity+'.'+electrolyte,
        throat_diffusivity=throat_diffusivity+'.'+electrolyte,
        conduit_lengths=conduit_lengths,
        conduit_shape_factors=conduit_shape_factors,
        pore_pressure=pore_pressure,
        pore_potential=pore_potential,
        throat_hydraulic_conductance=throat_hydraulic_conductance,
        throat_diffusive_conductance=(throat_diffusive_conductance + '.' +
                                      electrolyte),
        throat_valence=throat_valence+'.'+electrolyte,
        throat_temperature=throat_temperature,
        s_scheme=s_scheme)
