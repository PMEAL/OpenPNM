import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()


__all__ = [
    "water_correlation",
    "liquid_pure",
]


@docstr.dedent
def water_correlation(target, temperature='pore.temperature', salinity='pore.salinity'):
    r"""
    Calculates surface tension of pure water or seawater at atmospheric
    pressure using Eq. (28) given by Sharqawy et al. Values at
    temperature higher than the normal boiling temperature are calculated at
    the saturation pressure.

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    salinity : str
        The dictionary key containing the salinity values. Salinity must be
        expressed in g of salt per kg of solution (ppt).

    Returns
    -------
    value : ndarray
        A numpy ndarray containing surface tension of seawater in [N/m]

    Notes
    -----
    T must be in K, and S in g of salt per kg of phase, or ppt (parts per
    thousand). The correlation is valid for 273 < T < 313 K and
    0 < S < 40 g/kg within 0.2% accuracy.

    References
    ----------
    Sharqawy M. H., Lienhard J. H., and Zubair, S. M., Desalination and
    Water Treatment, 2010.

    """
    T = target[temperature]
    if salinity in target.keys():
        S = target[salinity]
    else:
        S = 0
    sigma_w = 0.2358*((1-(T/647.096))**1.256)*(1-0.625*(1-(T/647.096)))
    a1 = 2.2637334337E-04
    a2 = 9.4579521377E-03
    a3 = 3.3104954843E-02
    TC = T-273.15
    sigma_sw = sigma_w*(1+(a1*TC+a2)*np.log(1+a3*S))
    value = sigma_sw
    return value


@docstr.dedent
def liquid_pure(
    target,
    temperature='pore.temperature',
    critical_temperature='param.critical_temperature',
    boiling_temperature='param.boiling_temperature',
    critical_pressure='param.critical_pressure',
):
    r"""
    Uses Brock_Bird model

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.phase.T)s
    critical_temperature : str
        The dictionary key containing the critical temperature values (K)

    Returns
    -------
    value : ndarray
        A numpy ndarray containing surface tension values scaled to the
        temperature [N/m]

    """
    # brock_bird
    T = target[temperature]
    Tc = target[critical_temperature]
    Tr = T/Tc
    Tb = target[boiling_temperature]
    Tbr = Tb/Tc
    Pc = target[critical_pressure]/100000
    Q = 0.1196*(1 + Tbr*np.log(Pc/1.01325)/(1-Tbr))-0.279
    sigma = 1e-3 * Q * (Pc**(2/3)) * (Tc**(1/3)) * ((1-Tr)**(11/9))
    return sigma


def liquid_mixture(
    target,
    surface_tension='pore.surface_tension.*',
    density='pore.density.*',
    molecular_weight='param.molecular_weight.*o',
):
    r"""
    """
    sigmas = target.get_comp_vals(surface_tension)
    xs = target['pore.mole_fraction']
    rhos = target.get_comp_vals(density)  # kg/m3
    MWs = target.get_comp_vals(molecular_weight)  # g/mol
    rhoms = {k: rhos[k]/(MWs[k]/1000) for k in xs.keys()}  # mol/m3
    rhom_mix = np.vstack([rhoms[k]*xs[k] for k in xs.keys()]).sum(axis=0)
    sigma = 0.0
    for i, ki in enumerate(xs.keys()):
        for j, kj in enumerate(xs.keys()):
            num = xs[ki]*xs[kj]*(rhom_mix**2)*(sigmas[ki]*sigmas[kj])**0.5
            denom = rhoms[ki]*rhoms[kj]
            sigma += num/denom
    return sigma


if __name__ == "__main__":
    import chemicals as chem
    import openpnm as op
    from numpy.testing import assert_allclose

    pn = op.network.Demo()

    cbz = op.phase.Species(network=pn, species='chlorobenzene')
    cbz.add_model(propname='pore.surface_tension',
                  model=liquid_pure)
    cbz.add_model(propname='pore.density',
                  model=op.models.phase.density.liquid_pure)
    cbz['pore.molar_density'] = cbz['pore.density']/(cbz['param.molecular_weight']/1000)
    s_calc = cbz['pore.surface_tension'][0]
    s_ref = chem.interface.Brock_Bird(
        T=cbz['pore.temperature'][0],
        Tb=cbz['param.boiling_temperature'],
        Tc=cbz['param.critical_temperature'],
        Pc=cbz['param.critical_pressure'],
    )
    assert_allclose(s_ref, s_calc, rtol=1e-10, atol=0)

    bnz = op.phase.Species(network=pn, species='benzene')
    bnz.add_model(propname='pore.surface_tension',
                  model=liquid_pure)
    bnz.add_model(propname='pore.density',
                  model=op.models.phase.density.liquid_pure)
    bnz['pore.molar_density'] = bnz['pore.density']/(bnz['param.molecular_weight']/1000)
    s_calc = bnz['pore.surface_tension'][0]
    s_ref = chem.interface.Brock_Bird(
        T=bnz['pore.temperature'][0],
        Tb=bnz['param.boiling_temperature'],
        Tc=bnz['param.critical_temperature'],
        Pc=bnz['param.critical_pressure'],
    )
    assert_allclose(s_ref, s_calc, rtol=1e-10, atol=0)

    mix = op.phase.LiquidMixture(network=pn, components=[bnz, cbz])
    mix.x(bnz, 0.8)
    mix.x(cbz, 0.2)
    mix.add_model(propname='pore.surface_tension',
                  model=liquid_mixture)
    s_calc = mix['pore.surface_tension'][0]
    s_ref = chem.interface.Winterfeld_Scriven_Davis(
        xs=np.vstack(list(mix['pore.mole_fraction'].values()))[:, 0],
        sigmas=np.vstack(list(mix.get_comp_vals('pore.surface_tension').values()))[:, 0],
        rhoms=np.vstack(list(mix.get_comp_vals('pore.molar_density').values()))[:, 0],
    )
    assert_allclose(s_ref, s_calc, rtol=1e-2, atol=0)
