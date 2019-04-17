# from collections import ChainMap  # Might use eventually
import numpy as np
from openpnm.phases import GenericPhase as GenericPhase
from openpnm.utils import logging
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
        super().__init__(**kwargs)
        self.settings.update(settings)

        # Add any supplied phases to the phases list
        for comp in components:
            self.settings['components'].append(comp.name)
            self['pore.mole_fraction.'+comp.name] = np.nan

        self['pore.mole_fraction.all'] = np.zeros(self.Np, dtype=float)

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
        prop = '.'.join(key.split('.')[:2])
        invalid_keys = set(self.props()).difference(set(self.keys()))
        # invalid_keys = []
        # [invalid_keys.extend(item.keys()) for item in self.components.values()]
        if prop in invalid_keys:
            raise Exception(prop + ' already assigned to a component object')
        super().__setitem__(key, value)

    def props(self, incl_components=False, **kwargs):
        temp = []
        if incl_components:
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

    def _update_molfrac(self):
        # Update concentration.all
        self['pore.mole_fraction.all'] = 0.0
        dict_ = list(self['pore.mole_fraction'].values())
        if len(dict_) > 1:
            self['pore.mole_fraction.all'] = np.sum(dict_, axis=0)
        self['throat.mole_fraction.all'] = 0.0
        dict_ = list(self['throat.mole_fraction'].values())
        if len(dict_) > 1:
            self['throat.mole_fraction.all'] = np.sum(dict_, axis=0)

    def set_mole_fraction(self, component, values=[]):
        r"""
        Specify occupancy of a phase in each pore and/or throat

        Parameters
        ----------
        components : OpenPNM Phase object or name string
            The phase whose mole fraction is being specified

        values : array_like
            The mole fraction of ``component `` in each pore.  This array must
            be *Np*-long, with one value between 0 and 1 for each pore in the
            network.  If a scalar is received it is applied to all pores.

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
            logger.warning('Received Pvals contain mole fractions outside ' +
                           'the range of 0 to 1')
        if Pvals.size:
            self['pore.mole_fraction.' + component.name] = Pvals
        self._update_molfrac()

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
                self._update_molfrac()
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
