from openpnm.models.collections.phase import water
from openpnm.phase import Phase, _fetch_chemical_props
from openpnm.utils import Docorator


docstr = Docorator()


__all__ = [
    'Water',
]


@docstr.dedent
class Water(Phase):
    r"""
    Creates Phase object with preset values for Water

    Parameters
    ----------
    %(Phase.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        from thermo import Chemical
        a = Chemical('h2o')
        temp = _fetch_chemical_props(a)
        self.params.update(temp)
        self.models.update(water)
        self.regenerate_models()
