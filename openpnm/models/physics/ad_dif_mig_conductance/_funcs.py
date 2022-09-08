import numpy as _np
from openpnm.models import _doctxt


__all__ = ["ad_dif_mig"]


@_doctxt
def ad_dif_mig(phase,
               pore_pressure="pore.pressure",
               pore_potential="pore.potential",
               throat_hydraulic_conductance="throat.hydraulic_conductance",
               throat_diffusive_conductance="throat.diffusive_conductance",
               throat_valence="throat.valence",
               throat_temperature="throat.temperature",
               ion="",
               s_scheme="powerlaw"):
    r"""
    Calculate the advective-diffusive-migrative conductance of conduits
    in network

    Parameters
    ----------
    %(phase)s
    pore_pressure : str
        %(dict_blurb)s pore pressure
    pore_potential : str
        %(dict_blurb)s pore potential
    throat_hydraulic_conductance : str
        %(dict_blurb)s hydraulic conductance
    throat_diffusive_conductance : str
        %(dict_blurb)s diffusive conductance
    throat_valence : str
        %(dict_blurb)s throat ionic species valence
    pore_temperature : str
        %(dict_blurb)s pore temperature
    throat_temperature : str
        %(dict_blurb)s throat temperature
    ion : str
        Name of the ionic species
    s_scheme : str
        Name of the space discretization scheme to use

    Returns
    -------
    %(return_arr)s advection-diffusion-migration conductance

    Notes
    -----
    This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    'shape_factor' depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    # Add ion suffix to properties
    throat_diffusive_conductance = throat_diffusive_conductance + "." + ion
    throat_valence = throat_valence + "." + ion

    network = phase.project.network
    cn = network["throat.conns"]
    T = phase[throat_temperature]
    # Check if pressure and potential values exist, otherwise, assign zeros
    try:
        P = phase[pore_pressure]
    except KeyError:
        P = _np.zeros(shape=phase.Np, dtype=float)
    try:
        V = phase[pore_potential]
    except KeyError:
        V = _np.zeros(shape=phase.Np, dtype=float)
    z = phase[throat_valence]
    F = 96485.3329
    R = 8.3145
    gh = phase[throat_hydraulic_conductance]
    gd = phase[throat_diffusive_conductance]
    # .T below is for when gd is (Nt, 2) instead of (Nt, 1)
    gm = (gd.T * (z * F) / (R * T)).T
    delta_V = _np.diff(V[cn], axis=1).squeeze()
    delta_V = _np.append(delta_V, -delta_V)

    # Normal treatment when gd is Nt by 1
    if gd.size == network.Nt:
        gd = _np.tile(gd, 2)
        gm = _np.tile(gm, 2)
    # Special treatment when gd is not Nt by 1 (ex. mass partitioning)
    elif gd.size == 2 * network.Nt:
        gd = gd.reshape(network.Nt * 2, order="F")
        gm = gm.reshape(network.Nt * 2, order="F")
    else:
        raise Exception(f"Shape of {throat_diffusive_conductance} must either"
                        r" be (Nt,1) or (Nt,2)")

    # Migration
    mig = gm * delta_V

    # Advection
    Qij = -gh * _np.diff(P[cn], axis=1).squeeze()
    Qij = _np.append(Qij, -Qij)

    # Advection-migration
    adv_mig = Qij - mig

    # Peclet numbers
    Peij_adv_mig = adv_mig / gd  # includes advection and migration
    Peij_adv = Qij / gd  # includes advection only
    Peij_mig = mig / gd  # includes migration only
    # Filter values
    Peij_adv_mig[(Peij_adv_mig < 1e-10) & (Peij_adv_mig >= 0)] = 1e-10
    Peij_adv_mig[(Peij_adv_mig > -1e-10) & (Peij_adv_mig <= 0)] = -1e-10
    Peij_adv[(Peij_adv < 1e-10) & (Peij_adv >= 0)] = 1e-10
    Peij_adv[(Peij_adv > -1e-10) & (Peij_adv <= 0)] = -1e-10
    Peij_mig[(Peij_mig < 1e-10) & (Peij_mig >= 0)] = 1e-10
    Peij_mig[(Peij_mig > -1e-10) & (Peij_mig <= 0)] = -1e-10

    # Corrected advection-migration
    adv_mig = Peij_adv_mig * gd

    if s_scheme == "upwind":
        w = gd + _np.maximum(0, -adv_mig)
    elif s_scheme == "hybrid":
        w = _np.maximum(0, _np.maximum(-adv_mig, gd - adv_mig / 2))
    elif s_scheme == "powerlaw":
        w = gd * _np.maximum(
            0, (1 - 0.1 * _np.absolute(Peij_adv_mig)) ** 5
        ) + _np.maximum(0, -adv_mig)
    elif s_scheme == "powerlaw_upwind":
        w = (
            gd * _np.maximum(0, (1 - 0.1 * _np.absolute(Peij_adv)) ** 5)
            + _np.maximum(0, -Qij)
        ) + _np.maximum(0, mig)
    elif s_scheme == "exponential":
        w = -adv_mig / (1 - _np.exp(Peij_adv_mig))
    else:
        raise Exception("Unrecognized discretization scheme: " + s_scheme)
    w = w.reshape(network.Nt, 2, order="F")
    return w
