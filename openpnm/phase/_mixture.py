import numpy as np
import logging
from openpnm.phase import Phase
from openpnm.models.collections.phase import liquid_mixture
from openpnm.models.collections.phase import gas_mixture, binary_gas_mixture
from openpnm.models.phase.mixtures import mole_summation
from openpnm.utils import HealthDict
from openpnm.utils import Docorator, Workspace


logger = logging.getLogger(__name__)
docstr = Docorator()
ws = Workspace()


__all__ = [
    'Mixture',
    'BinaryGasMixture',
    'MultiGasMixture',
    'LiquidMixture',
]


class MixtureSettings:
    ...


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
            elif key.endswith('*'):
                key = key[:-1]
                if key.endswith('.'):
                    key = key[:-1]
                vals = self._get_comp_vals(key)
            else:
                raise KeyError(key)
        return vals

    def _get_comp_vals(self, propname):
        vals = {}
        for comp in self.components.keys():
            vals[comp] = self[propname + '.' + comp]
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

    def add_comp(self, component, mole_fraction=0.0):
        self['pore.mole_fraction.' + component.name] = mole_fraction

    def remove_comp(self, component):
        if hasattr(component, 'name'):
            component = component.name
        try:
            del self['pore.mole_fraction.' + component]
        except KeyError:
            pass

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


class LiquidMixture(Mixture):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.update(liquid_mixture())
        self.regenerate_models()

    def x(self, compname=None, x=None):
        r"""
        Helper method for getting and setting mole fractions of a component

        Parameters
        ----------
        compname : str, optional
            The name of the component, i.e. ``obj.name``.  If ``x`` is not
            provided this will *return* the mole fraction of the requested
            component.  If ``x`` is provided this will *set* the mole fraction
            of the specified component.
        x : scalar or ndarray, optional
            The mole fraction of the given species in the mixture. If not
            provided this method works as a *getter* and will return the
            mole fraction of the requested component.  If ``compname`` is not
            provided then the mole fractions of all components will be returned
            as a dictionary with the components names as keys.

        Notes
        -----
        This method is equivalent to
        ``mixture['pore.mole_fraction.<compname>'] = x``
        """
        if x is None:
            if compname is None:
                return self['pore.mole_fraction']
            else:
                self['pore.mole_fraction' + '.' + compname]
        else:
            self['pore.mole_fraction.' + compname] = x


class GasMixture(Mixture):

    def y(self, compname=None, y=None):
        r"""
        Helper method for getting and setting mole fractions of a component

        Parameters
        ----------
        compname : str, optional
            The name of the component i.e. ``obj.name``.  If ``y`` is not
            provided this will *return* the mole fraction of the requested
            component.  If ``y`` is provided this will *set* the mole fraction
            of the specified component.
        y : scalar or ndarray, optional
            The mole fraction of the given species in the mixture. If not
            provided this method works as a *getter* and will return the
            mole fraction of the requested component. If ``compname`` is also
            not provided then the mole fractions of all components will be
            returned as a dictionary with the components names as keys.

        Notes
        -----
        This method is equivalent to
        ``mixture['pore.mole_fraction.<compname>'] = y``
        """
        if y is None:
            if compname is None:
                return self['pore.mole_fraction']
            else:
                return self['pore.mole_fraction' + '.' + compname]
        else:
            self['pore.mole_fraction.' + compname] = y


class MultiGasMixture(GasMixture):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.clear()
        self.models.update(gas_mixture())
        self.regenerate_models()

class BinaryGasMixture(GasMixture):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.models.clear()
        self.models.update(binary_gas_mixture())
        self.regenerate_models()

    def add_comp(self, component, mole_fraction=0.0):
        if len(self['pore.mole_fraction'].keys()) >= 2:
            raise Exception("Binary mixtures cannot have more than 2 components"
                            + ", remove one first")
        super().add_component(component=component, mole_fraction=mole_fraction)


if __name__ == '__main__':
    import openpnm as op
    pn = op.network.Demo()
    o2 = op.phase.GasByName(network=pn, species='o2', name='pure_O2')
    n2 = op.phase.GasByName(network=pn, species='n2', name='pure_N2')
    air = op.phase.BinaryGasMixture(network=pn, components=[o2, n2], name='air')
    air.y(o2.name, 0.21)
    air['pore.mole_fraction.pure_N2'] = 0.79
    air.regenerate_models()
    print(air)

    ch4 = op.phase.GasByName(network=pn, species='ch4', name='methane')
    h2 = op.phase.GasByName(network=pn, species='h2', name='hydrogen')
    h2o = op.phase.GasByName(network=pn, species='h2o', name='water')
    co2 = op.phase.GasByName(network=pn, species='co2', name='co2')
    syngas = op.phase.MultiGasMixture(network=pn, components=[ch4, h2, h2o, co2],
                                      name='syngas')
    syngas.y(h2.name, 0.25)
    syngas.y(ch4.name, 0.25)
    syngas.y(h2o.name, 0.25)
    syngas.y(co2.name, 0.25)
    syngas.regenerate_models()
    print(syngas)































