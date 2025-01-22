import numpy as np
from pnmlib.core import get_data


__all__ = [
    "hydraulic_conductance",
    "hydraulic_conductance2",
]


def hydraulic_conductance(
    target,
    network,
    viscosity='viscosity',
    diameter='size',
):
    P1, P2 = network['throat.conns'][:].T
    pore_dia = network[f'pore.{diameter}'][:]
    throat_dia = network[f'throat.{diameter}'][:]
    pore_visc = target[f'pore.{viscosity}'][:]
    throat_visc = target[f'throat.{viscosity}'][:]
    gp1 = 8*np.pi*(pore_dia[P1]/2)**4/(8*pore_visc[P1])
    gp2 = 8*np.pi*(pore_dia[P2]/2)**4/(8*pore_visc[P2])
    gt = 8*np.pi*(throat_dia/2)**4/(8*throat_visc)
    gh = (1/gp1 + 1/gt + 1/gp2)**-1
    return gh


def hydraulic_conductance2(
    target,
    pore_viscosity,
    throat_viscosity,
    pore_diameter='network/pore.diameter',
    throat_diameter='network/throat.diameter',
):
    P1, P2 = get_data(target, 'network/throat.conns').T
    pore_dia = get_data(target, pore_diameter)
    throat_dia = get_data(target, throat_diameter)
    pore_visc = get_data(target, pore_viscosity)
    throat_visc = get_data(target, throat_viscosity)
    gp1 = 8*np.pi*(pore_dia[P1]/2)**4/(8*pore_visc[P1])
    gp2 = 8*np.pi*(pore_dia[P2]/2)**4/(8*pore_visc[P2])
    gt = 8*np.pi*(throat_dia/2)**4/(8*throat_visc)
    gh = (1/gp1 + 1/gt + 1/gp2)**-1
    return gh
