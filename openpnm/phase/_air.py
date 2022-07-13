from openpnm.models.collections.phase import air
from openpnm.phase import Phase, _fetch_chemical_props
from openpnm.utils import Docorator
from thermo import Mixture


docstr = Docorator()


__all__ = [
    'Air',
]


@docstr.dedent
class Air(Phase):
    r"""
    Creates a Phase object with preset models and values for air

    Parameters
    ----------
    %(Phase.parameters)s

    References
    ----------
    The correlations and constants for this class are taken from:

    ::

        E.W. Lemmon and R.T. Jacobsen, "Viscosity and Thermal Conductivity
        Equations for Nitrogen, Oxygen, Argon, and Air", Int. J. of
        Thermophysics, Vol. 25, No. 1, January 2004, pp. 21-69

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        a = Mixture(IDs=['o2', 'n2'], zs=[0.21, 0.79])
        temp = _fetch_chemical_props(a)
        self.params.update(temp)
        self.models.update(air(regen_mode='deferred'))
        # self.regenerate_models()
