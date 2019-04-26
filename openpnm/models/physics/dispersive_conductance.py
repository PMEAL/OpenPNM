r"""

.. autofunction:: openpnm.models.physics.dispersive_conductance.dispersion
.. autofunction:: openpnm.models.physics.ad_dif_conductance.generic_conductance

"""

from .ad_dif_conductance import generic_conductance


def dispersion(target,
               pore_area='pore.area',
               throat_area='throat.area',
               pore_diffusivity='pore.diffusivity',
               throat_diffusivity='throat.diffusivity',
               conduit_lengths='throat.conduit_lengths',
               conduit_shape_factors='throat.poisson_shape_factors',
               pore_pressure='pore.pressure',
               throat_hydraulic_conductance='throat.hydraulic_conductance',
               throat_diffusive_conductance='throat.diffusive_conductance',
               s_scheme='powerlaw'):
    r"""
    Calculate the dispersive conductance of conduits in network, where
    a conduit is ( 1/2 pore - full throat - 1/2 pore ). See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.

    pore_area : string
        Dictionary key of the pore area values

    throat_area : string
        Dictionary key of the throat area values

    pore_diffusivity : string
        Dictionary key of the pore diffusivity values

    throat_diffusivity : string
        Dictionary key of the throat diffusivity values

    conduit_lengths : string
        Dictionary key of the conduit length values

    conduit_shape_factors : string
        Dictionary key of the conduit DIFFUSION shape factor values

    pore_pressure : string
        Dictionary key of the pore pressure values

   throat_hydraulic_conductance : string
       Dictionary key of the throat hydraulic conductance values

   throat_diffusive_conductance : string
       Dictionary key of the throat diffusive conductance values

   s_scheme : string
       Name of the space discretization scheme to use

    Returns
    -------
    g : ndarray
        Array containing dispersive conductance values for conduits in the
        geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper flow_shape_factor argument.

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
