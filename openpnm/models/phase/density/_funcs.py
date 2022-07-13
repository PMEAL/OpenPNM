from openpnm.utils import Docorator
from chemicals import Vm_to_rho
import numpy as np


docstr = Docorator()


__all__ = [
    "ideal_gas",
    "water_correlation",
    "liquid_mixture",
    "liquid_pure",
    "mass_to_molar",
]


def ideal_gas(
    target,
    pressure='pore.pressure',
    temperature='pore.temperature',
    molecular_weight='param.molecular_weight',
):
    r"""
    Uses ideal gas law to calculate the mass density of an ideal gas

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    %(models.phase.P)s
    mol_weight : str
        Name of the dictionary key on ``target`` where the array containing
        molecular weight values is stored

    Returns
    -------
    %(models.phase.density.returns)s

    """
    P = target[pressure]
    T = target[temperature]
    try:
        # If target is a pure species, it should have molecular weight in params
        MW = target[molecular_weight]
    except KeyError:
        # Otherwise, get the mole weighted average value
        MW = target.get_mix_vals(molecular_weight)
    R = 8.31447  # J/(mol.K)
    value = P/(R*T)*MW/1000  # Convert to kg/m3
    return value


def water_correlation(
    target,
    temperature='pore.temperature',
    salinity='pore.salinity',
):
    r"""
    Calculates density of pure water or seawater at atmospheric pressure
    using Eq. (8) given by Sharqawy et. al [1]. Values at temperature higher
    than the normal boiling temperature are calculated at the saturation
    pressure.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    salinity : str
        Name of the dictionary key on ``target`` where the array containing
        salinity values is stored.

    Returns
    -------
    %(models.phase.density.returns)s

    Notes
    -----
     T must be in K, and S in g of salt per kg of phase, or ppt (parts per
        thousand)
    VALIDITY: 273 < T < 453 K; 0 < S < 160 g/kg;
    ACCURACY: 0.1 %

    References
    ----------
    [1] Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
    Water Treatment, 2010.

    """
    T = target[temperature]
    if salinity in target.keys():
        S = target[salinity]
    else:
        S = 0
    a1 = 9.9992293295E+02
    a2 = 2.0341179217E-02
    a3 = -6.1624591598E-03
    a4 = 2.2614664708E-05
    a5 = -4.6570659168E-08
    b1 = 8.0200240891E-01
    b2 = -2.0005183488E-03
    b3 = 1.6771024982E-05
    b4 = -3.0600536746E-08
    b5 = -1.6132224742E-11
    TC = T-273.15
    rho_w = a1 + a2*TC + a3*TC**2 + a4*TC**3 + a5*TC**4
    d_rho = b1*S + b2*S*TC + b3*S*(TC**2) + b4*S*(TC**3) + b5*(S**2)*(TC**2)
    rho_sw = rho_w + d_rho
    value = rho_sw
    return value


def liquid_mixture(
    target,
    temperature='pore.temperature',
    molecular_weight='param.molecular_weight',
    critical_temperature='param.critical_temperature',
    critical_volume='param.critical_volume',
    acentric_factor='param.acentric_factor',
):
    r"""
    Computes the density of a liquid mixture using the COrrospoding STAtes
    Liquid Density (COSTALD) method.

    Parameters
    ----------


    Returns
    -------
    density : ndarray
        The density of the liquid mixture in units of kg/m3.  Note that
        ``chemicals.volume.COSTALD`` returns molar volume, so this function
        converts it to density using the mole fraction weighted molecular
        weight of the mixture: :math:`MW_{mix} = \Sigma x_i \cdot MW_i`.

    Notes
    -----
    This is same approach used by ``chemicals.volume.COSTALD_mixture`` and
    exact numerical correspondance is confirmed. Unlike the ``chemicals``
    version, this function is vectorized over the conditions rather than
    the compositions, since a typical simulation has millions of pores each
    representing an independent conditions, but a mixture typically only has
    a few components.

    """
    # Fetch parameters for each pure component
    Tcs = target.get_comp_vals(critical_temperature)
    Vcs = target.get_comp_vals(critical_volume)
    omegas = target.get_comp_vals(acentric_factor)
    Xs = target['pore.mole_fraction']
    # Compute mixture values
    omegam = np.vstack([Xs[k]*omegas[k] for k in Xs.keys()]).sum(axis=0)
    Vm1 = np.vstack([Xs[k]*Vcs[k] for k in Xs.keys()]).sum(axis=0)
    Vm2 = np.vstack([Xs[k]*(Vcs[k])**(2/3) for k in Xs.keys()]).sum(axis=0)
    Vm3 = np.vstack([Xs[k]*(Vcs[k])**(1/3) for k in Xs.keys()]).sum(axis=0)
    Vm = 0.25*(Vm1 + 3*Vm2*Vm3)
    Tcm = 0.0
    for i, ki in enumerate(Xs.keys()):
        inner = 0.0
        for j, kj in enumerate(Xs.keys()):
            inner += Xs[ki]*Xs[kj]*(Vcs[ki]*Tcs[ki]*Vcs[kj]*Tcs[kj])**0.5
        Tcm += inner
    Tcm = Tcm/Vm
    # Convert molar volume to normal mass density
    MWs = target.get_comp_vals('param.molecular_weight')
    MWm = np.vstack([Xs[k]*MWs[k] for k in Xs.keys()]).sum(axis=0)
    T = target[temperature]
    rhoL = liquid_pure(
        target=target,
        temperature=T,
        molecular_weight=MWm,
        critical_temperature=Tcm,
        critical_volume=Vm,
        acentric_factor=omegam,
    )
    return rhoL


def liquid_pure(
    target,
    temperature='pore.temperature',
    critical_temperature='param.critical_temperature',
    critical_volume='param.critical_volume',
    acentric_factor='param.acentric_factor',
    molecular_weight='param.molecular_weight',
):
    r"""
    """
    Vc = target[critical_volume]
    Tc = target[critical_temperature]
    omega = target[acentric_factor]
    T = target[temperature]
    Tr = T/Tc
    V0 = 1 - 1.52816*(1-Tr)**(1/3) + 1.43907*(1-Tr)**(2/3) - 0.81446*(1-Tr) + \
        0.190454*(1-Tr)**(4/3)
    V1 = (-0.296123 + 0.386914*Tr - 0.0427258*Tr**2 - 0.0480645*Tr**3)/(Tr - 1.00001)
    Vs = Vc*V0*(1-omega*V1)
    MW = target[molecular_weight]
    rhoL = Vm_to_rho(Vm=Vs, MW=MW)
    return rhoL


def mass_to_molar(
    target,
    molecular_weight='param.molecular_weight',
    density='pore.density',
):
    r"""
    Calculates the molar density from the molecular weight and mass density

    Parameters
    ----------
    %(models.target.parameters)s
    mol_weight : str
        The dictionary key containing the molecular weight in kg/mol
    density : str
        The dictionary key containing the density in kg/m3

    Returns
    -------
    value : ndarray
        A numpy ndrray containing molar density values [mol/m3]

    """
    MW = target[molecular_weight]
    rho = target[density]
    value = rho/MW
    return value


if __name__ == "__main__":

    import chemicals as chem
    import openpnm as op
    from numpy.testing import assert_allclose

    pn = op.network.Demo()

    h2o = op.phase.Species(network=pn, species='water')
    h2o.add_model(propname='pore.density',
                  model=op.models.phase.density.liquid_pure)
    Vm = chem.COSTALD(
        T=h2o['pore.temperature'][0],
        Tc=h2o['param.critical_temperature'],
        Vc=h2o['param.critical_volume'],
        omega=h2o['param.acentric_factor'],
    )
    rho_ref = Vm_to_rho(Vm, h2o['param.molecular_weight'])
    rho_calc = h2o['pore.density'][0]
    assert_allclose(rho_ref, rho_calc, rtol=1e-10, atol=0)

    etoh = op.phase.Species(network=pn, species='ethanol')
    etoh.add_model(propname='pore.density',
                   model=op.models.phase.density.liquid_pure)

    vodka = op.phase.LiquidMixture(network=pn, components=[h2o, etoh])
    vodka.x(h2o.name, 0.5)
    vodka.x(etoh.name, 0.5)
    vodka.add_model(propname='pore.density',
                    model=op.models.phase.density.liquid_mixture)
    Vm = chem.COSTALD_mixture(
        T=vodka['pore.temperature'][0],
        xs=np.vstack(list(vodka['pore.mole_fraction'].values()))[:, 0],
        Tcs=list(vodka.get_comp_vals('param.critical_temperature').values()),
        Vcs=list(vodka.get_comp_vals('param.critical_volume').values()),
        omegas=list(vodka.get_comp_vals('param.acentric_factor').values()),
    )
    rho_ref = Vm_to_rho(Vm, vodka.get_mix_vals('param.molecular_weight')[0])
    rho_calc = vodka['pore.density'][0]
    assert_allclose(rho_ref, rho_calc, rtol=1e-10, atol=0)























