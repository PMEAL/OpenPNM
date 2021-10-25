r"""
Pore-scale models for calculating the advective-diffusive-migrative
conductance of conduits.
"""
import numpy as _np

__all__ = ["ad_dif_mig"]


def ad_dif_mig(
    target,
    pore_pressure="pore.pressure",
    pore_potential="pore.potential",
    throat_hydraulic_conductance="throat.hydraulic_conductance",
    throat_diffusive_conductance="throat.diffusive_conductance",
    throat_valence="throat.valence",
    pore_temperature="pore.temperature",
    throat_temperature="throat.temperature",
    ion="",
    s_scheme="powerlaw",
):
    r"""
    Calculate the advective-diffusive-migrative conductance of conduits
    in network, where a conduit is ( 1/2 pore - full throat - 1/2 pore ).
    See the notes section.

    Parameters
    ----------
    target : OpenPNM Object
        The object which this model is associated with. This controls the
        length of the calculated array, and also provides access to other
        necessary properties.
    pore_pressure : string
        Dictionary key of the pore pressure values
    pore_potential : string
        Dictionary key of the pore potential values
    throat_hydraulic_conductance : string
        Dictionary key of the throat hydraulic conductance values
    throat_diffusive_conductance : string
        Dictionary key of the throat diffusive conductance values
    throat_valence : string
        Dictionary key of the throat ionic species valence values
    pore_temperature : string
        Dictionary key of the pore temperature values
    throat_temperature : string
        Dictionary key of the throat temperature values
    ion : string
        Name of the ionic species
    s_scheme : string
        Name of the space discretization scheme to use

    Returns
    -------
    g : ndarray
        Array containing advective-diffusive-migrative conductance values for
        conduits in the geometry attached to the given physics object.

    Notes
    -----
    (1) This function requires that all the necessary phase properties already
    be calculated.

    (2) This function calculates the specified property for the *entire*
    network then extracts the values for the appropriate throats at the end.

    (3) This function assumes cylindrical throats with constant cross-section
    area. Corrections for different shapes and variable cross-section area can
    be imposed by passing the proper conduit_shape_factors argument.

    (4) shape_factor depends on the physics of the problem, i.e. diffusion-like
    processes and fluid flow need different shape factors.

    """
    # Add ion suffix to properties
    throat_diffusive_conductance = throat_diffusive_conductance + "." + ion
    throat_valence = throat_valence + "." + ion

    network = target.project.network
    throats = network.map_throats(throats=target.Ts, origin=target)
    phase = target.project.find_phase(target)
    cn = network["throat.conns"][throats]
    T = phase[throat_temperature][throats]
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
    if gd.size == throats.size:
        gd = _np.tile(gd, 2)
        gm = _np.tile(gm, 2)
    # Special treatment when gd is not Nt by 1 (ex. mass partitioning)
    elif gd.size == 2 * throats.size:
        gd = gd.reshape(throats.size * 2, order="F")
        gm = gm.reshape(throats.size * 2, order="F")
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
    w = w.reshape(throats.size, 2, order="F")
    return w
