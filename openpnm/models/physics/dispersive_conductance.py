r"""

.. autofunction:: openpnm.models.physics.diffusive_conductance.ordinary_diffusion

"""

from .misc import generic_conductance


def dispersion(target,
               pore_area='pore.area',
               throat_area='throat.area',
               pore_diffusivity='pore.diffusivity',
               pore_pressure='pore.pressure',
               throat_hydraulic_conductance='throat.hydraulic_conductance',
               throat_diffusive_conductance='throat.diffusive_conductance',
               throat_diffusivity='throat.diffusivity',
               conduit_lengths='throat.conduit_lengths',
               conduit_shape_factors='throat.poisson_shape_factors',
               s_scheme='powerlaw'):
    r"""

    """
    return generic_conductance(
        target=target,
        transport_type='dispersion',
        pore_area=pore_area,
        throat_area=throat_area,
        pore_diffusivity=pore_diffusivity,
        throat_diffusivity=throat_diffusivity,
        conduit_lengths=conduit_lengths,
        conduit_shape_factors=conduit_shape_factors,
        pore_pressure=pore_pressure,
        throat_hydraulic_conductance=throat_hydraulic_conductance,
        throat_diffusive_conductance=throat_diffusive_conductance,
        s_scheme=s_scheme)
