from numpy import pi as _pi
from openpnm.utils import Docorator


docstr = Docorator()
docstr.params['models.geometry.throat_cross_sectional_area.returns'] = \
    r"""areas : ndarray
            A numpy ndarray containing throat cross-sectional area values"""


@docstr.dedent
def cylinder(target, throat_diameter='throat.diameter'):
    r"""
    Calculate throat cross-sectional area for a cylindrical throat

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tdia.parameters)s

    Returns
    -------
    %(models.geometry.throat_cross_sectional_area.returns)s

    """
    diams = target[throat_diameter]
    value = _pi / 4 * diams**2
    return value


@docstr.dedent
def cuboid(target, throat_diameter='throat.diameter'):
    r"""
    Calculate throat cross-sectional area for a cuboid throat

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tdia.parameters)s

    Returns
    -------
    %(models.geometry.throat_cross_sectional_area.returns)s

    """
    diams = target[throat_diameter]
    value = (diams)**2
    return value


def rectangle(target, throat_diameter='throat.diameter'):
    r"""
    Calculate throat cross-sectional area for a rectangular throat

    Parameters
    ----------
    %(models.target.parameters)s
    %(models.geometry.tdia.parameters)s

    Returns
    -------
    %(models.geometry.throat_cross_sectional_area.returns)s

    """
    return target[throat_diameter]
