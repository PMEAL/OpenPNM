import numpy as np
from openpnm.phase import GenericPhase
import openpnm.models.phase as mods
from openpnm.utils import HealthDict, PrintableList, SubDict
from openpnm.utils import Docorator, SettingsAttr
from openpnm.utils import logging


logger = logging.getLogger(__name__)
docstr = Docorator()


class GenericMixtureSettings:
    r"""
    The following settings are specific to Mixture objects

    Parameters
    ----------
    components : list of strings
        The names of each pure component object that constitute the mixture

    """

    components = []
    prefix = 'mix'


@docstr.get_sections(base='GenericMixture', sections=['Parameters'])
@docstr.dedent
class GenericMixture(GenericPhase):
    r"""
    Creates Phase object that represents a multicomponent mixture system
    consisting of a given list of GenericPhases as components.

    Parameters
    ----------
    %(GenericPhase.parameters)s
    components : list
        A list of all components that constitute this mixture

    Notes
    -----
    All mixtures assume that mole fractions are always stored as
    ``'pore.mole_fraction'`` and that concentration is always stored as
    ``'pore.concentration'``.

    """

    def __init__(self, components=[], settings=None, **kwargs):
        self.settings = SettingsAttr(GenericMixtureSettings, settings)
        super().__init__(settings=self.settings, **kwargs)

        # Add any supplied phases to the phases list
        for comp in components:
            self.set_component(comp)
            # self.set_mole_fraction(comp, np.nan)

        self.add_model(propname='pore.mole_fraction.all',
                       model=mods.mixtures.mole_summation)

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
        lines = super().__str__()
        lines = '\n'.join((lines, 'Component Phases', horizontal_rule))
        for item in self.components.values():
            lines = '\n'.join((lines, item.__module__.replace('__', '')
                               + ' : ' + item.name))
        lines = '\n'.join((lines, horizontal_rule))
        return lines

    def _del_comps(self, components):
        if not isinstance(components, list):
            components = [components]
        # Remove from settings:
        for comp in components:
            self.settings['components'].remove(comp.name)
        # Remove data from dict:
        for item in list(self.keys()):
            if item.endswith(comp.name):
                self.pop(item)
        self['pore.mole_fraction.all'] = np.nan

    def _get_comps(self):
        comps = {}
        comps.update({item: self.project[item]
                      for item in sorted(self.settings['components'])})
        return comps

    def _set_comps(self, components):
        if not isinstance(components, list):
            components = [components]
        temp = self.settings['components'].copy()
        temp.extend([val.name for val in components])
        comps = list(set(temp))
        self.settings['components'] = comps
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


class LiquidMixture(GenericMixture):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_model(propname='pore.molecular_weight',
                       model=mods.mixture.mixture_molecular_weight,)
        self.add_model(propname='pore.viscosity',
                       model=mods.viscosity.liquid_mixture_viscosity,)
        self.add_model(propname='pore.critical_volume',
                       model=mods.critical_properties.liquid_mixture_critical_volume,)
        self.add_model(propname='pore.critical_temperature',
                       model=mods.critical_properties.liquid_mixture_critical_temperature,)
        self.add_model(propname='pore.acentric_factor',
                       model=mods.critical_properties.mixture_acentric_factor,)
        self.add_model(propname='pore.density',
                       model=mods.density.liquid_mixture_density,)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.thermal_conductivity.liquid_mixture_thermal_conductivity,)
        self.add_model(propname='pore.heat_capacity',
                       model=mods.heat_capacity.mixture_heat_capacity,)


class GasMixture(GenericMixture):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_model(propname='pore.molecular_weight',
                       model=mods.mixtures.mixture_molecular_weight,)
        self.add_model(propname='pore.gas_viscosity',
                       model=mods.viscosity.gas_mixture_viscosity,)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.thermal_conductivity.gas_mixture_thermal_conductivity,)
        self.add_model(propname='pore.heat_capacity',
                       model=mods.heat_capacity.mixture_heat_capacity,)
        self.add_model(propname='pore.LJ_epsilon',
                       model=mods.diffusivity.gas_mixture_LJ_epsilon,)
        self.add_model(propname='pore.LJ_sigma',
                       model=mods.diffusivity.gas_mixture_LJ_sigma,)
        self.add_model(propname='pore.LJ_omega',
                       model=mods.diffusivity.gas_mixture_LJ_collision_integral,)
        self.add_model(propname='pore.diffusivity',
                       model=mods.diffusivity.gas_mixture_diffusivity,)
