import numpy as _np
from openpnm.models import _doctxt


__all__ = ["ad_dif"]


@_doctxt
def ad_dif(
    phase,
    pore_pressure='pore.pressure',
    throat_hydraulic_conductance='throat.hydraulic_conductance',
    throat_diffusive_conductance='throat.diffusive_conductance',
    s_scheme='powerlaw'
):
    r"""
    Calculates the advective-diffusive conductance of conduits in network.

    Parameters
    ----------
    %(phase)s
    pore_pressure : str
        %(dict_blurb)s pore pressure
    throat_hydraulic_conductance : str
        %(dict_blurb)s hydraulic conductance
    throat_diffusive_conductance : str
        %(dict_blurb)s throat diffusive conductance
    s_scheme : str
        Name of the space discretization scheme to use

    Returns
    -------
    %(return_arr)s advection-diffuvsion conductance values

    Notes
    -----
    This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the
    end.

    This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area
    can be imposed by passing the proper conduit_shape_factors argument
    when computing the diffusive and hydraulic conductances.

    shape_factor depends on the physics of the problem, i.e.
    diffusion-like processes and fluid flow need different shape factors.

    """
    network = phase.project.network
    cn = network['throat.conns']
    # Find g for half of pore 1, throat, and half of pore 2
    P = phase[pore_pressure]
    gh = phase[throat_hydraulic_conductance]
    gd = phase[throat_diffusive_conductance]
    if gd.size == network.Nt:
        gd = _np.tile(gd, 2)
    # Special treatment when gd is not Nt by 1 (ex. mass partitioning)
    elif gd.size == 2 * network.Nt:
        gd = gd.reshape(network.Nt * 2, order='F')
    else:
        raise Exception(f"Shape of {throat_diffusive_conductance} must either"
                        r" be (Nt,1) or (Nt,2)")

    Qij = -gh * _np.diff(P[cn], axis=1).squeeze()
    Qij = _np.append(Qij, -Qij)

    Peij = Qij / gd
    Peij[(Peij < 1e-10) & (Peij >= 0)] = 1e-10
    Peij[(Peij > -1e-10) & (Peij <= 0)] = -1e-10

    # Correct the flow rate
    Qij = Peij * gd

    if s_scheme == 'upwind':
        w = gd + _np.maximum(0, -Qij)
    elif s_scheme == 'hybrid':
        w = _np.maximum(0, _np.maximum(-Qij, gd - Qij / 2))
    elif s_scheme == 'powerlaw':
        w = gd * _np.maximum(0, (1 - 0.1 * _np.absolute(Peij))**5) + \
            _np.maximum(0, -Qij)
    elif s_scheme == 'exponential':
        w = -Qij / (1 - _np.exp(Peij))
    else:
        raise Exception('Unrecognized discretization scheme: ' + s_scheme)
    w = w.reshape(network.Nt, 2, order='F')
    return w
