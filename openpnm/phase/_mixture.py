import numpy as np
import logging
from copy import deepcopy
from openpnm.phase import Phase
from openpnm.models.collections.phase import liquid_mixture, gas_mixture
from openpnm.models.phase.mixtures import mole_summation
from openpnm.utils import HealthDict, PrintableList, SubDict
from openpnm.utils import Docorator, Workspace


logger = logging.getLogger(__name__)
docstr = Docorator()
ws = Workspace()


__all__ = [
    'Mixture',
    'GasMixture',
    'LiquidMixture',
]


class MixtureSettings:
    r"""
    The following settings are specific to Mixture objects

    Parameters
    ----------
    components : list of strings
        The names of each pure component object that constitute the mixture

    """
    components = []


@docstr.get_sections(base='Mixture', sections=['Parameters'])
@docstr.dedent
class Mixture(Phase):
    r"""
    Creates Phase object that represents a multicomponent mixture
    consisting of a given list of individual phase objects as components.

    Parameters
    ----------
    %(Phase.parameters)s
    components : list
        A list of all components that constitute this mixture

    """

    def __init__(self, components=[], name='mixture_?', **kwargs):
        self._components = []
        super().__init__(name=name, **kwargs)
        self.settings._update(MixtureSettings())

        # Add any supplied phases to the list components
        self.components = components

    def __getitem__(self, key):
        try:
            vals = super().__getitem__(key)
        except KeyError:
            # If key ends in component name, fetch it
            if key.split('.')[-1] in self.components.keys():
                comp = self.project[key.split('.')[-1]]
                vals = comp[key.rsplit('.', maxsplit=1)[0]]
            else:
                raise KeyError(key)
        return vals

    def __str__(self):
        horizontal_rule = '―' * 78
        temp = super().__str__().split(horizontal_rule)
        lines = ''
        lines = '\n'.join((lines, 'Component Phases', horizontal_rule))
        for item in self.components.values():
            module = self.__module__
            module = ".".join([x for x in module.split(".") if not x.startswith("_")])
            cname = module + '.' + item.__class__.__name__
            lines = '\n'.join((lines, cname + ' : ' + item.name))
        temp.insert(2, lines + '\n')
        lines = horizontal_rule.join(temp)
        return lines

    def _get_comps(self):
        comps = {}
        comps.update({item: self.project[item]
                      for item in sorted(self['pore.mole_fraction'].keys())})
        return comps

    def _set_comps(self, components):
        if not isinstance(components, list):
            components = [components]
        # Add mole_fraction array to dict, filled with nans
        for item in components:
            if 'pore.mole_fraction.' + item.name not in self.keys():
                self['pore.mole_fraction.' + item.name] = 0.0

    components = property(fget=_get_comps, fset=_set_comps)

    def set_x(self, compname, x):
        r"""
        Helper method for setting mole fraction of a component

        Parameters
        ----------
        compname : str
            The name of the component object, i.e. ``obj.name``
        x : scalar or ndarray
            The mole fraction of the given species in the mixture.

        Notes
        -----
        This method is equivalent to
        ``mixture['pore.mole_fraction.<compname>'] = x``
        """
        self['pore.mole_fraction.' + compname] = x

    def get_x(self, compname):
        r"""
        Helper method for retrieving the mole fraction of a component

        Parameters
        ----------
        compname : str
            The name of the component object, i.e. ``obj.name``

        Notes
        -----
        This method is equivalent to ``mixture['pore.mole_fraction.<compname>']``
        """
        return self['pore.mole_fraction.' + compname]

    def check_mixture_health(self):
        r"""
        Checks the "health" of the mixture

        Calculates the mole fraction of all species in each pore and returns
        an list of where values are too low or too high

        Returns
        -------
        health : dict
            A HealthDict object containing lists of locations where the mole
            fractions are not unity. One value indicates locations that are
            too high, and another where they are too low.

        """
        h = HealthDict()
        h['mole_fraction_too_low'] = []
        h['mole_fraction_too_high'] = []
        conc = mole_summation(target=self)
        lo = np.where(conc < 1.0)[0]
        hi = np.where(conc > 1.0)[0]
        if len(lo) > 0:
            h['mole_fraction_too_low'] = lo
        if len(hi) > 0:
            h['mole_fraction_too_high'] = hi
        return h


class ComponentHandler:

    def __iadd__(self, component):
        target = self._find_target()
        if 'pore.mole_fraction.' + component.name not in target.keys():
            target['pore.mole_fraction.'+component.name] = 0.0

    def __isub__(self, component):
        target = self._find_target()
        del target['pore.mole_fraction.'+component.name]

    def _find_target(self):
        """
        Finds and returns the parent object to self.
        """
        for proj in ws.values():
            for obj in proj:
                if hasattr(obj, "components"):
                    if obj.components is self:
                        return obj


class LiquidMixture(Mixture):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(liquid_mixture())
        self.regenerate_models()


class GasMixture(Mixture):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(gas_mixture())
        self.regenerate_models()


if __name__ == '__main__':
    import openpnm as op
    pn = op.network.Demo()
    o2 = op.phase.GasByName(network=pn, species='o2', name='pure_O2')
    n2 = op.phase.GasByName(network=pn, species='n2', name='pure_N2')
    air = op.phase.GasMixture(network=pn, components=[o2, n2])
    # air.set_x(o2.name, 0.21)
    # air['pore.mole_fraction.pure_N2'] = 0.79
    air.regenerate_models()
    print(air)
