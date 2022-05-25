from openpnm.models.collections.phase import air
from openpnm.phase import GenericPhase
from openpnm.utils import Docorator


docstr = Docorator()


@docstr.dedent
class Air(GenericPhase):
    r"""
    Creates a Phase object with preset models and values for air

    Parameters
    ----------
    %(GenericPhase.parameters)s

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
        self.models.update(air())
        self.regenerate_models()
