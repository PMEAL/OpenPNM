from numpy import pi


def cylinders_in_series(target,
                        pore_diameter='pore.diameter',
                        throat_diameter='throat.diameter',
                        n_cylinders=5,
                        throat_length=None,
                        return_elements=False):
    r"""

    """
    network = target.network


def spheres_and_cylinders(target,
                         pore_diameter='pore.diameter',
                         throat_diameter='throat.diameter',
                         throat_length=None,
                         return_elements=False):
    r"""

    """
    network = target.network


def pyramids_and_cuboids(target,
                        pore_diameter='pore.diameter',
                        throat_diameter='throat.diameter',
                        throat_length=None,
                        return_elements=False):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns]/2).T
    Rt = target[throat_diameter]/2


def cones_and_cylinders(target,
                       pore_diameter='pore.diameter',
                       throat_diameter='throat.diameter',
                       throat_length='throat.length',
                       return_elements=False):
    r"""

    """
    network = target.network
    R1, R2 = (target[pore_diameter][network.conns]/2).T
    Rt = target[throat_diameter]/2
    Lt = target[throat_length]
    L1 = Lt - R1
    L2 = Lt - R2

    alpha1 = (R1 - Rt)/L1
    beta1 = 1 / (1/(Rt**3) - 1/(R1**3))
    alpha2 = (R2-Rt)/L2
    beta2 = 1 / (1/(Rt**3) - 1/(R2**3))
    g1 = (3*alpha1*pi/8) * beta1
    g2 = (3*alpha2*pi/8) * beta2
    gt = pi*Rt**4/(8*Lt)

    if return_elements:
        g = {'pore1': g1, 'throat': gt, 'pore2': g2}
    else:
        g = (1/g1 + 1/gt + 1/g2)**-1
    return g
