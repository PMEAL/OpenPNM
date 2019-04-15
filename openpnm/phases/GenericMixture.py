from openpnm.phases import GenericPhase as GenericPhase
import openpnm.models.phases as _models
import numpy as np
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
            self['throat.mole_fraction.'+comp.name] = np.nan

        self['pore.mole_fraction.all'] = np.zeros(self.Np, dtype=float)
        self['throat.mole_fraction.all'] = np.zeros(self.Nt, dtype=float)

        self.pop('pore.temperature', None)
        self.pop('pore.pressure', None)

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

    def set_mole_fraction(self, component, Pvals=[], Tvals=[]):
        r"""
        Specify occupancy of a phase in each pore and/or throat

        Parameters
        ----------
        components : OpenPNM Phase object or name string
            The phase whose mole fraction is being specified

        Pvals : array_like
            The mole fraction of ``component `` in each pore.  This array must
            be *Np*-long, with one value between 0 and 1 for each pore in the
            network.  If a scalar is received it is applied to all pores.

        Tvals : array_like
            The mole fraction of ``component`` in each throat.  This array must
            be *Nt*-long, with one value between 0 and 1, for each throat in
            the network.  If a scalar is received it is applied to all throats.

        """
        if type(component) == str:
            component = self.components[component]
        Pvals = np.array(Pvals, ndmin=1)
        Tvals = np.array(Tvals, ndmin=1)
        if component not in self.project:
            raise Exception(f"{component.name} doesn't belong to this project")
        else:
            if component.name not in self.settings['components']:
                self.settings['components'].append(component.name)
        if np.any(Pvals > 1.0) or np.any(Pvals < 0.0):
            logger.warning('Received Pvals contain mole fractions outside ' +
                           'the range of 0 to 1')
        if np.any(Tvals > 1.0) or np.any(Tvals < 0.0):
            logger.warning('Received Tvals contain mole fractions outside ' +
                           'the range of 0 to 1')
        if Pvals.size:
            self['pore.mole_fraction.' + component.name] = Pvals
        if Tvals.size:
            self['throat.mole_fraction.' + component.name] = Tvals
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
        if np.any(self[element + '.mole_fraction.all'] != 1.0):
            self._update_molfrac()
            if np.any(self[element + '.mole_fraction.all'] != 1.0):
                raise Exception('Mole fraction does not add to unity in all ' +
                                element + 's')
        vals = np.zeros([self._count(element=element)], dtype=float)
        try:
            for comp in self.components.values():
                vals += comp[prop]*self[element + '.mole_fraction.' + comp.name]
        except KeyError:
            vals = super().interleave_data(prop)
        return vals
