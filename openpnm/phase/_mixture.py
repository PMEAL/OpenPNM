import logging

import numpy as np

from openpnm.models.collections.phase import (
    binary_gas_mixture,
    standard_gas_mixture,
    standard_liquid_mixture,
)
from openpnm.models.phase.mixtures import mole_summation
from openpnm.phase import Phase
from openpnm.utils import (
    Docorator,
    HealthDict,
    Workspace,
    get_printable_labels,
    get_printable_props,
)

logger = logging.getLogger(__name__)
docstr = Docorator()
ws = Workspace()


__all__ = [
    'Mixture',
    'BinaryGas',
    'GasMixture',
    'LiquidMixture',
    'StandardGasMixture',
    'StandardLiquidMixture',
]


class MixtureSettings:
    ...


@docstr.get_sections(base='Mixture', sections=['Parameters'])
@docstr.dedent
class Mixture(Phase):
    r"""
    Creates ``Phase`` object that represents a multicomponent mixture
    consisting of a given list of individual ``Species`` objects as
    components.

    Parameters
    ----------
    %(Phase.parameters)s
    components : list of ``Species`` objects
        A list of all components that constitute the mixture

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
            # If key ends in component name, fetch it directly
            if key.split('.')[-1] in self.components.keys():
                comp = self.project[key.split('.')[-1]]
                vals = comp[key.rsplit('.', maxsplit=1)[0]]
            elif key[-1] in ['*', '.']:  # Fetch values for ALL components
                key = key[:-1]
                while key.endswith('.'):  # If ended with .* also remove .
                    key = key[:-1]
                vals = self.get_comp_vals(key)
            else:
                raise KeyError(key)
        return vals

    def regenerate_models(self, **kwargs):
        for c in self.components.values():
            c.regenerate_models(**kwargs)
        super().regenerate_models(**kwargs)

    def get_comp_vals(self, propname):
        r"""
        Get a dictionary of the requested values from each component in the
        mixture

        Parameters
        ----------
        propname : str
            The property to fetch, such as ``'pore.viscosity'``.

        Returns
        -------
        vals : dict
            A dictionary with each component name as the key and the requested
            property as the value.

        """
        if not isinstance(propname, str):
            return propname
        if propname.endswith('*'):
            try:
                return self[propname]
            except KeyError:
                pass
        try:
            vals = {}
            for comp in self.components.values():
                vals[comp.name] = comp[propname]
            return vals
        except KeyError:
            msg = f'{propname} not found on at least one component'
            raise Exception(msg)

    def get_mix_vals(self, propname, mode='linear', power=1):
        r"""
        Get the mole fraction weighted value of a given property for the mixture

        Parameters
        ----------
        propname : str
            The property to fetch, such as ``'pore.viscosity'``.
        mode : str
            The type of weighting to use. Options are:

            ============== ====================================================
            mode
            ============== ====================================================
            'linear'       (default) Basic mole fraction weighting of the form:
                           :math:`z = \Sigma (x_i \cdot \z_i)`
            'logarithmic'  Uses the natural logarithm of the property as:
                           :math:`ln(z) = \Sigma (x_i \cdot ln(\z_i))`
            'power         Applies an exponent to the property as:
                           :math:`\z^{power} = \Sigma (x_i \cdot \z_i^{power})`
            ============== ====================================================

        power : scalar
            If ``mode='power'`` this indicates the value of the exponent,
            otherwise this is ignored.

        Returns
        -------
        vals : ndarray
            An ndarray with the requested values obtained using the specified
            weigthing.

        """
        Xs = self['pore.mole_fraction']
        ys = self.get_comp_vals(propname)
        if mode == 'linear':
            z = np.vstack([Xs[k]*ys[k] for k in Xs.keys()]).sum(axis=0)
        elif mode == 'logarithmic':
            z = np.vstack([Xs[k]*np.log(ys[k]) for k in Xs.keys()]).sum(axis=0)
            z = np.exp(z)
        elif mode == 'power':
            z = np.vstack([Xs[k]*(ys[k])**power for k in Xs.keys()]).sum(axis=0)
            z = z**(1/power)
        return z

    @property
    def info(self):  # pragma: no cover
        hr = '―' * 78
        lines = '\n' + 'Component Phases'
        for item in self.components.values():
            lines += '\n' + "═"*78 + '\n' + item.__repr__() + '\n' + hr + '\n'
            lines += get_printable_props(item, suffix=item.name)
            lines += '\n' + hr + '\n'
            lines += get_printable_labels(item, suffix=item.name)
        temp = super().__str__().split(hr)
        temp.insert(-1, lines + '\n')
        lines = hr.join(temp)
        print(lines)

    def _get_comps(self):
        comps = {}
        for k in self.keys():
            if k.startswith('pore.mole_fraction'):
                name = k.split('.')[-1]
                comps[name] = self.project[name]
        return comps

    def _set_comps(self, components):
        if not isinstance(components, list):
            components = [components]
        # Add mole_fraction array to dict, filled with nans
        for item in components:
            if 'pore.mole_fraction.' + item.name not in self.keys():
                self['pore.mole_fraction.' + item.name] = np.nan

    components = property(fget=_get_comps, fset=_set_comps)

    def add_comp(self, component, mole_fraction=np.nan):
        r"""
        Helper method to add a component

        Parameters
        ----------
        component : Species object
            The pure ``Species`` object to add
        mole_fraction : float or array_like
            The mole fraction of the given component in each pore. Can either
            be a scalar which is applied to every pore, or an Np-long
            ``ndarray`` with values for each pore.

        Notes
        -----
        This method just adds `'pore.mole_fraction.<component.name>'`
        to the `Mixture` dictionary, and assigns the given `mole_fraction`.
        """
        self['pore.mole_fraction.' + component.name] = mole_fraction

    def remove_comp(self, component):
        r"""
        Helper method to remove a component from the mixture

        Parameters
        ----------
        component : Species object or str name
            The `Species` to remove from the mixture. Can either be a handle to
            the object, or the object's `name`.

        Notes
        -----
        This method just calls
        `del mixture['pore.mole_fraction.<component.name>']` to remove the
        given component.

        """
        if hasattr(component, 'name'):
            component = component.name
        try:
            del self['pore.mole_fraction.' + component]
        except KeyError:
            pass

    def check_mixture_health(self):
        r"""
        Checks the state of health of the mixture

        Calculates the mole fraction of all species in each pore and returns
        an list of where values are too low or too high,

        Returns
        -------
        health : dict
            A HealthDict object containing lists of locations where the mole
            fractions are not unity. One value indicates locations that are
            too high, and another where they are too low. This `dict` evaluates
            to `True` if all items are healthy.

        """
        h = HealthDict()
        h['mole_fraction_too_low'] = []
        h['mole_fraction_too_high'] = []
        conc = mole_summation(phase=self)
        lo = np.where(conc < 1.0)[0]
        hi = np.where(conc > 1.0)[0]
        if len(lo) > 0:
            h['mole_fraction_too_low'] = lo
        if len(hi) > 0:
            h['mole_fraction_too_high'] = hi
        return h


class LiquidMixture(Mixture):
    r"""
    Wrapper class for creating a liquid mixture.

    This class has a helper method `x` which allows for setting and getting
    of compositions.
    """

    def x(self, compname=None, x=None):
        r"""
        Helper method for getting and setting mole fractions of a component

        Parameters
        ----------
        compname : str, optional
            The name of the component, i.e. ``obj.name``.  If ``x`` is not
            provided this will *return* the mole fraction of the requested
            component.  If ``x`` is provided this will *set* the mole fraction
            of the specified component to ``x``.
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
        if hasattr(compname, 'name'):
            compname = compname.name
        if x is None:
            if compname is None:
                return self['pore.mole_fraction']
            else:
                self['pore.mole_fraction' + '.' + compname]
        else:
            self['pore.mole_fraction.' + compname] = x


class GasMixture(Mixture):
    r"""
    Wrapper class for creating a gas mixture.

    This class has a helper method `y` which allows for setting and getting
    of compositions.
    """

    def y(self, compname=None, y=None):
        r"""
        Helper method for getting and setting mole fractions of a component

        Parameters
        ----------
        compname : str, optional
            The name of the component i.e. ``obj.name``.  If ``y`` is not
            provided this will *return* the mole fraction of the requested
            component.  If ``y`` is provided this will *set* the mole fraction
            of the specified component to ``y``.
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
        if hasattr(compname, 'name'):
            compname = compname.name
        if y is None:
            if compname is None:
                return self['pore.mole_fraction']
            else:
                return self['pore.mole_fraction' + '.' + compname]
        else:
            self['pore.mole_fraction.' + compname] = y


class BinaryGas(GasMixture):
    r"""
    A ``GasMixture`` class that has a predefined set of models for computing the
    properties of the mixture as a function of composition.

    The number of components is capped at 2 (i.e. binary) since the calculation
    of diffusion coefficients using the Lennard-Jones method in the `chemicals`
    package only works for binary mixtures.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_model_collection(binary_gas_mixture)
        # Running models doesn't make sense as compositions are not set yet
        # self.regenerate_models()

    def __setitem__(self, key, value):
        if key.startswith('pore.mole_fraction'):
            if key not in self.keys():
                try:
                    d = self['pore.mole_fraction'].keys()
                    if len(d) >= 2:
                        msg = "Binary mixtures cannot have more than" \
                            " 2 components, remove one first"
                        raise Exception(msg)
                except KeyError:
                    pass
        super().__setitem__(key, value)


class StandardLiquidMixture(LiquidMixture):
    r"""
    A ``LiquidMixture`` class that also has a predefined suite of models for
    computing the properties of the mixture as a function of composition.

    The models are taken from the `chemicals` package and can be viewed with
    ``print(liq_mix.models)``, where `liq_mix` is an instance of this class.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_model_collection(standard_liquid_mixture)
        # Running models doesn't make sense as compositions are not set yet
        # self.regenerate_models()


class StandardGasMixture(GasMixture):
    r"""
    A ``GasMixture`` class that also has a predefined suite of models for
    computing the properties of the mixture as a function of composition.

    The models are taken from the `chemicals` package and can be viewed with
    ``print(gas_mix.models)``, where `gas_mix` is an instance of this class.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_model_collection(standard_gas_mixture)
        # Running models doesn't make sense as compositions are not set yet
        # self.regenerate_models()
