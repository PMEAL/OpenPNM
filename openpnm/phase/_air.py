from openpnm.models.collections.phase import air
from openpnm.phase import Phase, _fetch_chemical_props
from openpnm.phase import StandardGas, StandardGasMixture
from openpnm.utils import Docorator


docstr = Docorator()


__all__ = [
    'Air',
    '_AirMixture',
]


@docstr.dedent
class Air(Phase):
    r"""
    Creates a Phase object with preset models and values for air

    Parameters
    ----------
    %(Phase.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        from thermo import Mixture
        a = Mixture(IDs=['o2', 'n2'], zs=[0.21, 0.79])
        temp = _fetch_chemical_props(a)
        self.params.update(temp)
        self.models.update(air)
        self.regenerate_models()


class _AirMixture(StandardGasMixture):  # pragma: no cover

    def __init__(self, network, **kwargs):
        o2 = StandardGas(network=network, species='o2')
        n2 = StandardGas(network=network, species='n2')
        super().__init__(network=network, components=[o2, n2], **kwargs)
        self.y(o2, 0.21)
        self.y(n2, 0.79)
        self.regenerate_models()
