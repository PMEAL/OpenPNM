import numpy as np
from openpnm.utils import Docorator


docstr = Docorator()


__all__ = [
    "water_correlation",
    "liquid_pure_brock_bird",
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
def liquid_pure_brock_bird(
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
    T = target[temperature]
    Tc = target[critical_temperature]
    Tr = T/Tc
    Tb = target[boiling_temperature]
    Tbr = Tb/Tc
    Pc = target[critical_pressure]/100000
    Q = 0.1196*(1 + Tbr*np.log(Pc/1.01325)/(1-Tbr))-0.279
    sigma = 1e-3 * Q * (Pc**(2/3)) * (Tc**(1/3)) * ((1-Tr)**(11/9))
    return sigma


if __name__ == "__main__":
    import chemicals as chem
    import openpnm as op
    from numpy.testing import assert_allclose

    pn = op.network.Demo()

    cbz = op.phase.Species(network=pn, species='chlorobenzen')
    cbz.add_model(propname='pore.surface_tension',
                  model=liquid_pure_brock_bird)
    k_calc = cbz['pore.surface_tension'][0]
    k_ref = chem.interface.Brock_Bird(
        T=cbz['pore.temperature'][0],
        Tb=cbz['param.boiling_temperature'],
        Tc=cbz['param.critical_temperature'],
        Pc=cbz['param.critical_pressure'],
    )
    assert_allclose(k_ref, k_calc, rtol=1e-10, atol=0)


