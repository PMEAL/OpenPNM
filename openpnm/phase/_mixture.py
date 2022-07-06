import numpy as np
import logging
from openpnm.phase import Phase
from openpnm.models.collections.phase import liquid_mixture, gas_mixture
from openpnm.utils import HealthDict, PrintableList, SubDict
from openpnm.utils import Docorator


logger = logging.getLogger(__name__)
docstr = Docorator()


__all__ = [
    'Mixture',

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

    def __init__(self, components=[], name='mixture_#', **kwargs):
        self._components = []
        super().__init__(name=name, **kwargs)
        self.settings._update(MixtureSettings())

        # Add any supplied phases to the components
        for comp in components:
            self.set_component(comp)

    def __getitem__(self, key):
        try:
            vals = super().__getitem__(key)
        except KeyError:
            try:
                # If key ends in component name, fetch it
                if key.split('.')[-1] in self.settings['components']:
                    comp = self.project[key.split('.')[-1]]
                    vals = comp[key.rsplit('.', maxsplit=1)[0]]
                    return vals
                # If does not end in component name, see if components have it
                vals = SubDict()
                for comp in self.components.keys():
                    vals[key + '.' + comp] = self.components[comp][key]
            except KeyError:
                vals = self.interleave_data(key)
        return vals

    def __setitem__(self, key, value):
        # Prevent writing 'element.property.component' on mixture
        invalid_keys = set(self.props(deep=True)).difference(set(self.props()))
        if key in invalid_keys:
            raise Exception(key + ' already assigned to a component object')
        # If writing a concentration, use set_concentration setter
        if key.startswith('pore.concentration'):
            self.set_concentration(component=key.split('.')[-1], values=value)
        super().__setitem__(key, value)

    def props(self, deep=False, **kwargs):
        temp = PrintableList()
        if deep:
            for item in self.components.values():
                temp.extend([prop+'.'+item.name for prop in item.props(**kwargs)])
        temp.extend(super().props(**kwargs))
        temp.sort()
        return temp

    def __str__(self):
        horizontal_rule = '―' * 78
        temp = super().__str__().split(horizontal_rule)
        lines = ''
        lines = '\n'.join((lines, 'Component Phases', horizontal_rule))
        for item in self.components.values():
            lines = '\n'.join((lines, item.__module__.replace('__', '')
                               + ' : ' + item.name))
        temp.insert(2, lines + '\n')
        lines = horizontal_rule.join(temp)
        return lines

    def _del_comps(self, components):
        if not isinstance(components, list):
            components = [components]
        # Remove from settings:
        for comp in components:
            self._components.remove(comp.name)
        # Remove data from dict:
        for item in list(self.keys()):
            if item.endswith(comp.name):
                self.pop(item)
        self['pore.mole_fraction.all'] = np.nan

    def _get_comps(self):
        comps = {}
        comps.update({item: self.project[item]
                      for item in sorted(self._components)})
        return comps

    def _set_comps(self, components):
        if not isinstance(components, list):
            components = [components]
        temp = self._components.copy()
        temp.extend([val.name for val in components])
        comps = list(set(temp))
        self._components = comps
        # Add mole_fraction array to dict, filled with nans
        for item in comps:
            self['pore.mole_fraction.' + item] = np.nan

    components = property(fget=_get_comps, fset=_set_comps)

    def set_component(self, component, mode='add'):
        r"""
        Add another component to the mixture

        Parameters
        ----------
        component : GenericPhase
            The phase object of the component to add.  Can also be a list
            of phases to add more than one at a time.
        mode : str {'add', 'remove'}
            Indicates whether to add or remove the give component(s)

        """
        if mode == 'add':
            self._set_comps(component)
        if mode == 'remove':
            self._del_comps(component)

    def check_mixture_health(self):
        r"""
        Checks the "health" of the mixture

        Calculates the mole fraction of all species in each pore and returns
        an list of where values are too low or too high

        Returns
        -------
        health : dict
            A HealtDict object containing lists of locations where the mole
            fractions are not unity.  One value indicates locations that are
            too high, and another where they are too low.

        """
        h = HealthDict()
        h['mole_fraction_too_low'] = []
        h['mole_fraction_too_high'] = []
        self.regenerate_models('pore.mole_fraction.all')
        lo = np.where(self['pore.mole_fraction.all'] < 1.0)[0]
        hi = np.where(self['pore.mole_fraction.all'] > 1.0)[0]
        if len(lo) > 0:
            h['mole_fraction_too_low'] = lo
        if len(hi) > 0:
            h['mole_fraction_too_high'] = hi
        return h


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
