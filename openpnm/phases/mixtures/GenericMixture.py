# from collections import ChainMap  # Might use eventually
import numpy as np
from openpnm.phases import GenericPhase as GenericPhase
from openpnm.utils import logging, HealthDict, PrintableList
logger = logging.getLogger(__name__)


class GenericMixture(GenericPhase):
    r"""
    Creates Phase object that represents a multicomponent mixture system
    consisting of a given list of OpenPNM Phase objects as components.

    Parameters
    ----------
    network : OpenPNM Network object
        The network to which this phase object will be attached.

    components : list of OpenPNM Phase objects
        A list of all components that constitute this mixture

    project : OpenPNM Project object, optional
        The Project with which this phase should be associted.  If a
        ``network`` is given then this is ignored and the Network's project
        is used.  If a ``network`` is not given then this is mandatory.

    name : string, optional
        The name of the phase.  This is useful to keep track of the objects
        throughout the simulation.  The name must be unique to the project.
        If no name is given, one is generated.

    """
    def __init__(self, components=[], settings={}, **kwargs):
        self.settings.update({'components': [],
                              })
        super().__init__(settings={'prefix': 'mix'}, **kwargs)
        self.settings.update(settings)

        # Add any supplied phases to the phases list
        for comp in components:
            self.settings['components'].append(comp.name)
            self['pore.mole_fraction.'+comp.name] = 0.0

        self['pore.mole_fraction.all'] = np.zeros(self.Np, dtype=float)
        logger.warning('Mixtures are a beta feature and functionality may ' +
                       'change in future versions')

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
            lines = '\n'.join((lines, item.__module__.replace('__', '') +
                               ' : ' + item.name))
        lines = '\n'.join((lines, horizontal_rule))
        return lines

    def _update_total_molfrac(self):
        # Update mole_fraction.all
        self['pore.mole_fraction.all'] = 0.0
        dict_ = list(self['pore.mole_fraction'].values())
        if len(dict_) > 1:
            self['pore.mole_fraction.all'] = np.sum(dict_, axis=0)
        self['throat.mole_fraction.all'] = 0.0
        dict_ = list(self['throat.mole_fraction'].values())
        if len(dict_) > 1:
            self['throat.mole_fraction.all'] = np.sum(dict_, axis=0)

    def update_concentrations(self, mole_fraction='pore.mole_fraction'):
        r"""
        Re-calculates the concentration of each species in the mixture based
        on the current mole fractions.

        This method looks up the mole fractions *and* the density of the
        mixture, then finds the respective concentrations in $mol/m^{3}$.

        Parameters
        ----------

        """
        density = self['pore.molar_density']
        for item in self.components.values():
            mf = self['pore.mole_fraction.'+item.name]
            self['pore.concentration.'+item.name] = density*mf

    def update_mole_fractions(self, concentration='pore.concentration',
                              molar_density='pore.molar_density'):
        r"""
        Re-calculates mole fraction of each species in mixture based on the
        current concentrations.

        This method looks up the concentration of each species (using the
        optionally specified concentration dictionary key), and calculates
        the mole fraction.  Optionally, it can use a molar density for the
        mixture and N-1 concentrations to determine the Nth concentration and
        all species mole fractions.

        Parameters
        ----------
        concentration : string, optional
            The dictionary key pointing to the desired concentration values.
            The default is 'pore.concentration'.  Given this value, lookups
            are performed for each species in the mixture.
        molar_density : string, optional
            The dictionary key pointing to the molar density of the mixture.
            If not given (default), then 'pore.molar_density' is used. If
            there are N-1 concentrations specified, then ``molar_density`` is
            automatically used to find the Nth concentration.

        Notes
        -----
        The method does not return any values.  Instead it updates the mole
        fraction arrays of each species directly.
        """
        concentrations = [concentration + '.' + comp for comp
                          in self.settings['components']
                          if concentration + '.' + comp in self.keys()]
        if len(concentrations) == len(self.components):
            # Find total number of moles per unit volume
            density = 0.0
            for conc in concentrations:
                density += self[conc]
            # Normalize moles per unit volume for each species by the total
            for conc in concentrations:
                element, quantity, component = conc.split('.')
                self[element+'.mole_fraction.'+component] = self[conc]/density
        elif len(concentrations) == (len(self.components) - 1):
            # Find mole fraction of N-1 species
            mol_frac = 0.0
            density = self[molar_density]
            for conc in concentrations:
                element, quantity, component = conc.split('.')
                self[element+'.mole_fraction.'+component] = self[conc]/density
                mol_frac += self[element+'.mole_fraction.'+component]
            # Find mole fraction of Nth species using molar_density
            given_comps = [conc.split('.')[2] for conc in concentrations]
            all_comps = self.settings['components']
            component = list(set(all_comps).difference(set(given_comps)))[0]
            self[element+'.mole_fraction.'+component] = 1 - mol_frac
            # [self[element+'.concentration.'+component] = (1 - mol_frac)*density
        else:
            raise Exception('Insufficient mole fraction values found for ' +
                            'component species')

    def set_concentration(self, component, values=[]):
        r"""
        Specify mole fraction of each component in each pore

        Parameters
        ----------
        components : OpenPNM Phase object or name string
            The phase whose mole fraction is being specified

        values : array_like
            The concentration of the given ``component `` in each pore.  This
            array must be *Np*-long, with one value for each pore in the
            network.  If a scalar is received it is applied to all pores.

        See Also
        --------
        set_mole_fraction

        """
        if type(component) == str:
            component = self.components[component]
        Pvals = np.array(values, ndmin=1)
        if component not in self.project:
            raise Exception(f"{component.name} doesn't belong to this project")
        else:
            if component.name not in self.settings['components']:
                self.settings['components'].append(component.name)
        if np.any(Pvals < 0.0):
            logger.warning('Received values contain negative concentrations')
        if Pvals.size:
            self['pore.concentration.' + component.name] = Pvals

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

        See Also
        --------
        set_concentration

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
            logger.warning('Received values contain mole fractions outside ' +
                           'the range of 0 to 1')
        if Pvals.size:
            self['pore.mole_fraction.' + component.name] = Pvals
        self._update_total_molfrac()

    def _get_comps(self):
        comps = {item: self.project[item] for item in self.settings['components']}
        return comps

    def _set_comps(self, components):
        if not isinstance(components, list):
            components = [components]
        self.settings['components'] = [val.name for val in components]

    components = property(fget=_get_comps, fset=_set_comps)

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
                    raise Exception('Mole fraction does not add to unity in all ' +
                                    element + 's')
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
