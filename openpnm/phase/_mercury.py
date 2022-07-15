from openpnm.models.collections.phase import mercury
from openpnm.phase import Phase, _fetch_chemical_props
from thermo import Chemical


__all__ = [
    'Mercury',
]


class Mercury(Phase):
    r"""
    Creates Phase object with and preset values and pore-scale models for
    mercury

    Parameters
    ----------
    %(Phase.parameters)s

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        a = Chemical('hg')
        temp = _fetch_chemical_props(a)
        self.params.update(temp)
        self.add_model_collection(mercury())
        self.regenerate_models()
