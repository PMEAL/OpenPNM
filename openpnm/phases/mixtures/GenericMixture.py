# from collections import ChainMap  # Might use eventually
import numpy as np
from openpnm.phases import GenericPhase as GenericPhase
from openpnm.utils import HealthDict, PrintableList, Docorator, GenericSettings
from openpnm.utils import logging
logger = logging.getLogger(__name__)
docstr = Docorator()


class GenericMixtureSettings(GenericSettings):
    r"""
    The following settings are specific to Mixture objects

    Parameters
    ----------
    components : list of strings
        The names of each pure component object that constitute the mixture

    """

    components = []
    prefix = 'mix'


@docstr.get_sectionsf('GenericMixture', sections=['Parameters'])
@docstr.dedent
class GenericMixture(GenericPhase):
    r"""
    Creates Phase object that represents a multicomponent mixture system
    consisting of a given list of OpenPNM Phase objects as components.

    Parameters
    ----------
    %(GenericPhase.parameters)s
    components : list of OpenPNM Phase objects
        A list of all components that constitute this mixture

    Notes
    -----
    All mixtures assume that mole fractions are always stored as
    ``'pore.mole_fraction'`` and that concentration is always stored as
    ``'pore.concentration'``.

    """

    def __init__(self, components=[], settings={}, **kwargs):
        self.settings._update_settings_and_docs(GenericMixtureSettings)
        self.settings.update(settings)
        super().__init__(settings={'prefix': 'mix'}, **kwargs)

        # Add any supplied phases to the phases list
        for comp in components:
            self.set_component(comp)
            self.set_mole_fraction(comp, np.nan)
        self['pore.mole_fraction.all'] = np.nan

        logger.warning('Mixtures are a beta feature and functionality may '
                       + 'change in future versions')

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
                else:
                    raise KeyError
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

    def _update_total_molfrac(self):
        # Update mole_fraction.all
        self.pop('pore.mole_fraction.all', None)
        temp = np.zeros((self.Np, ), dtype=float)
        for item in self.keys():
            if item.startswith('pore.mole_fraction'):
                temp += self[item]
        self['pore.mole_fraction.all'] = temp

    def _reset_molfracs(self):
        for item in self.settings['components']:
            self['pore.mole_fraction.'+item] = np.nan
        self['pore.mole_fraction.all'] = np.nan

    def update_mole_fractions(self, free_comp=None):
        r"""
        Updates mole fraction values so the sum of all mole fractions is
        1.0 in each pore.

        Parameters
        ----------
        free_comp : OpenPNM object
            Specifies which component's mole fraction should be adjusted to
            enforce the sum of all mole fractions to equal 1.0.  If not given
            then a search is conducted to find a component whose mole
            fractions are *nans* and this component is treated as
            ``free_comp``.  If more than one component is found during
            this search an error is raised.  If *all* the components are
            *nans* then an attempt is made to determine mole fractions from
            the concentrations.

        """
        # Scan the list of components and find if any have nans
        comps = self.settings['components'].copy()
        hasnans = []
        for item in comps:
            if 'pore.mole_fraction.'+item not in self.keys():
                self['pore.mole_fraction.'+item] = np.nan
            if np.any(np.isnan(self['pore.mole_fraction.'+item])):
                hasnans.append(item)
        # If only one has nans then update it's mole fraction to make sum = 1.0
        if len(hasnans) == 1:
            comp = hasnans[0]
            self['pore.mole_fraction.'+comp] = 1.0
            comps.remove(comp)
            for item in comps:
                self['pore.mole_fraction.'+comp] -= self['pore.mole_fraction.'+item]
            self._update_total_molfrac()
        # If none have nans then either adjust given free_comp or raise error
        elif len(hasnans) == 0:
            if free_comp:
                self['pore.mole_fraction.'+free_comp] = np.nan
                self.update_mole_fractions()
            else:
                raise Exception('No free component found, specify which to adjust')
        # If they all have nans, then compute them all based on concentrations
        elif len(hasnans) == len(comps):
            total_conc = np.zeros((self.Np, ), dtype=float)
            for item in self.settings['components']:
                try:
                    total_conc += self['pore.concentration.'+item]
                except KeyError:
                    logger.warn('Concentrations for all species are not yet '
                                + ' defined')
            for item in self.settings['components']:
                mf = self['pore.concentration.'+item]/total_conc
                self['pore.mole_fraction.'+item] = mf
            self._update_total_molfrac()

    def set_concentration(self, component, values=[]):
        r"""
        Specify the concentration of a component in each pore

        Parameters
        ----------
        component : OpenPNM Phase object or name
            The phase object of the component whose concentration is being
            specified
        values : array_like
            The concentration of the given ``component `` in each pore.  This
            array must be *Np*-long, with one value for each pore in the
            network.  If a scalar is received it is applied to all pores.

        See Also
        --------
        set_mole_fraction

        Notes
        -----
        This method automatically sets all the mole fractions to *nans*.

        """
        if type(component) == str:
            component = self.components[component]
        Pvals = np.array(values, ndmin=1)
        if component not in self.components.values():
            raise Exception('Given component not part of mixture')
        if Pvals.size:
            key = 'pore.concentration.' + component.name
            super().__setitem__(key, Pvals)
        self._reset_molfracs()

    def set_mole_fraction(self, component, values=[]):
        r"""
        Specify mole fraction of each component in each pore

        Parameters
        ----------
        components : OpenPNM Phase object or name string
            The phase whose mole fraction is being specified
        values : array_like
            The mole fraction of the given ``component `` in each pore.  This
            array must be *Np*-long, with one value between 0 and 1 for each
            pore in the network.  If a scalar is received it is applied to
            all pores.

        """
        if type(component) == str:
            component = self.components[component]
        Pvals = np.array(values, ndmin=1)
        if component not in self.project:
            raise Exception(f"{component.name} doesn't belong to this project")
        else:
            if component.name not in self.settings['components']:
                self.settings['components'].append(component.name)
        if np.any(Pvals > 1.0) or np.any(Pvals < 0.0):
            logger.warning('Received values contain mole fractions outside '
                           + 'the range of 0 to 1')
        if Pvals.size:
            self['pore.mole_fraction.' + component.name] = Pvals
        self._update_total_molfrac()

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
        comps = {item: self.project[item] for item in self.settings['components']}
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
        component : OpenPNM Phase object
            The phase object of the component to add.  Can also be a list
            of phases to add more than one at a time.
        mode : string {'add', 'remove'}
            Indicates whether to add or remove the give component(s)

        """
        if mode == 'add':
            self._set_comps(component)
        if mode == 'remove':
            self._del_comps(component)

    def interleave_data(self, prop):
        r"""
        Gathers property values from component phases to build a single array

        If the requested ``prop`` is not on this Mixture, then a search is
        conducted on all associated components objects, and values from each
        are assembled into a single array.

        Parameters
        ----------
        prop : string
            The property to be retrieved

        Returns
        -------
        array : ND-array
            An array containing the specified property retrieved from each
            component phase and assembled based on the specified mixing rule

        """
        element = prop.split('.')[0]
        if element == 'pore':
            if np.any(self[element + '.mole_fraction.all'] != 1.0):
                self._update_total_molfrac()
                if np.any(self[element + '.mole_fraction.all'] != 1.0):
                    raise Exception('Mole fraction does not add to unity in all '
                                    + element + 's')
        vals = np.zeros([self._count(element=element)], dtype=float)
        try:
            for comp in self.components.values():
                vals += comp[prop]*self[element+'.mole_fraction.'+comp.name]
        except KeyError:
            vals = super().interleave_data(prop)
        return vals

    def check_mixture_health(self):
        r"""
        Checks the "health" of the mixture

        Calculates the mole fraction of all species in each pore and returns
        an list of where values are too low or too high

        Returns
        -------
        health : dict
            A HealtDict object containing lists of locations where the mole
            fractions are not unity.  One value indiates locations that are
            too high, and another where they are too low.

        """
        h = HealthDict()
        h['mole_fraction_too_low'] = []
        h['mole_fraction_too_high'] = []
        self._update_total_molfrac()
        lo = np.where(self['pore.mole_fraction.all'] < 1.0)[0]
        hi = np.where(self['pore.mole_fraction.all'] > 1.0)[0]
        if len(lo) > 0:
            h['mole_fraction_too_low'] = lo
        if len(hi) > 0:
            h['mole_fraction_too_high'] = hi
        return h
