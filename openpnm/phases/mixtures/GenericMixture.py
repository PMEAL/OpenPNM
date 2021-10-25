# from collections import ChainMap  # Might use eventually
import numpy as np
from chemicals import Vm_to_rho, rho_to_Vm
from chemicals import numba_vectorized
from chemicals import R, k
from openpnm.phases import GenericPhase as GenericPhase
from openpnm.utils import HealthDict, PrintableList, SubDict
from openpnm.utils import Docorator, GenericSettings
from openpnm import models
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


@docstr.get_sections(base='GenericMixture', sections=['Parameters'])
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

    def _update_total_molfrac(self):
        temp = np.zeros((self.Np, ), dtype=float)
        comps = list(self.components.values())
        for item in comps:
            temp += self['pore.mole_fraction.' + item.name]
        self['pore.mole_fraction.all'] = temp

    def _reset_molfracs(self):
        for item in self.settings['components']:
            self['pore.mole_fraction.'+item] = np.nan
        self['pore.mole_fraction.all'] = np.nan

    def update_concentrations(self, mode='all'):
        r"""
        Updates all unspecified concentrations from mole fractions and molar
        density.

        Parameters
        ----------
        mode : str {'all' (default), 'update'}
            If 'all' then all concentrations will be updated based on the
            current mole fractions.  If 'update', then only concentrations
            that are missing or filled with *nans* will be updated.

        Notes
        -----
        This function requires that molar density be defined on the mixture.

        """
        if 'pore.molar_density' not in self.keys():
            raise Exception("Cannot update concentration without "
                            + "\'pore.molar_density\'")
        comps = self.settings['components'].copy()
        for item in comps:
            c = self['pore.mole_fraction.' + item] * self['pore.molar_density']
            if mode == 'all':
                self['pore.concentration.' + item] = c
            elif mode == 'update':
                if 'pore.concentration.' + item not in self.items():
                    self['pore.concentration.' + item] = c

    def update_mole_fractions(self, free_comp=None):
        r"""
        Updates mole fraction values so the sum of all mole fractions is
        1.0 in each pore.

        Parameters
        ----------
        free_comp : OpenPNM object
            Specifies which component's mole fraction should be adjusted to
            enforce the sum of all mole fractions to equal 1.0.  See ``Notes``
            for more details.

        Notes
        -----
        If ``free_comp`` is not given then a search is conducted to find a
        component whose mole fractions are *nans* and this component is
        treated as ``free_comp``.

        If more than one component, or *no* components, is/are found during
        this search, an attempt is made to determine mole fractions from
        the concentrations, such that $x_i = C_i/C_total$.

        """
        # Parse the free_comp argument
        if hasattr(free_comp, 'name'):  # Get name if give openpnm object
            free_comp = free_comp.name
        if free_comp is not None:  # Set the component's mole fraction to nan
            self['pore.mole_fraction.' + free_comp] = np.nan
        # Scan the list of components and find if any have nans
        comps = self.settings['components'].copy()
        hasnans = []
        for item in comps:
            # If component not found, give it nans
            if 'pore.mole_fraction.'+item not in self.keys():
                self['pore.mole_fraction.'+item] = np.nan
            hasnans.append(np.any(np.isnan(self['pore.mole_fraction.' + item])))
        # If only one has nans then update it's mole fraction to make sum = 1.0
        if hasnans.count(True) == 1:
            # Remove free_comp from list
            comp = comps.pop(hasnans.index(True))
            # Set comp with nans to 1.0, then deduct other mole fracs
            self['pore.mole_fraction.'+comp] = 1.0
            for item in comps:
                self['pore.mole_fraction.' + comp] -= \
                    self['pore.mole_fraction.' + item]
        # If all or none of the components have values, then use concentrations
        else:
            # First find total concentration of all components
            total_conc = 0.0
            try:
                for item in comps:
                    total_conc += self['pore.concentration.' + item]
                # Then find mole fractions as C_i/C_total
                for item in comps:
                    mf = self['pore.concentration.' + item] / total_conc
                    self['pore.mole_fraction.' + item] = mf
            except KeyError:
                logger.warning("Insufficient concentrations defined, cannot "
                               + "recalculate mole fractions")
        self._update_total_molfrac()
        self.regenerate_models(['pore.molar_mass'])

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

        """
        if isinstance(component, str):
            component = self.components[component]
        Pvals = np.array(values, ndmin=1)
        if component not in self.components.values():
            raise Exception('Given component not part of mixture')
        if Pvals.size:
            key = 'pore.concentration.' + component.name
            super().__setitem__(key, Pvals)

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
        # Parse the input args
        if isinstance(component, str):
            component = self.components[component]
        Pvals = np.array(values, ndmin=1)
        # If given component not part of mixture, set it
        if component.name not in self.settings['components']:
            self.set_component(component)
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


class LiquidMixture(GenericMixture):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_model(propname='pore.molecular_weight',
                       model=mixture_molecular_weight,
                       regen_mode='deferred')
        self.add_model(propname='pore.viscosity',
                       model=liquid_mixture_viscosity,
                       regen_mode='deferred')
        self.add_model(propname='pore.critical_volume',
                       model=liquid_mixture_critical_volume,
                       regen_mode='deferred')
        self.add_model(propname='pore.critical_temperature',
                       model=liquid_mixture_critical_temperature,
                       regen_mode='deferred')
        self.add_model(propname='pore.acentric_factor',
                       model=mixture_acentric_factor,
                       regen_mode='deferred')
        self.add_model(propname='pore.density',
                       model=liquid_mixture_density,
                       regen_mode='deferred')
        self.add_model(propname='pore.thermal_conductivity',
                       model=liquid_mixture_thermal_conductivity,
                       regen_mode='deferred')
        self.add_model(propname='pore.heat_capacity',
                       model=mixture_heat_capacity,
                       regen_mode='deferred')


class GasMixture(GenericMixture):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.add_model(propname='pore.molecular_weight',
                       model=mixture_molecular_weight,
                       regen_mode='deferred')
        self.add_model(propname='pore.gas_viscosity',
                       model=gas_mixture_viscosity,
                       regen_mode='deferred')
        self.add_model(propname='pore.thermal_conductivity',
                       model=gas_mixture_thermal_conductivity,
                       regen_mode='deferred')
        self.add_model(propname='pore.heat_capacity',
                       model=mixture_heat_capacity,
                       regen_mode='deferred')
        self.add_model(propname='pore.LJ_epsilon',
                       model=gas_mixture_LJ_epsilon,
                       regen_mode='deferred')
        self.add_model(propname='pore.LJ_sigma',
                       model=gas_mixture_LJ_sigma,
                       regen_mode='deferred')
        self.add_model(propname='pore.LJ_omega',
                       model=gas_mixture_LJ_collision_integral,
                       regen_mode='deferred')
        self.add_model(propname='pore.diffusion_coefficient',
                       model=gas_mixture_diffusivity,
                       regen_mode='deferred')


def gas_mixture_viscosity(target):
    ys = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    MWs = [c.parameters['molecular_weight'] for c in target.components.values()]
    mus = [c['pore.viscosity'] for c in target.components.values()]
    # Should be a one-liner, but not working
    # mu = numba_vectorized.Herning_Zipperer(ys, mus, MWs, sqrtMWs)
    num = 0.0
    denom = 0.0
    for i in range(len(ys)):
        num += ys[i]*mus[i]*np.sqrt(MWs[i])
        denom += ys[i]*np.sqrt(MWs[i])
    mu = num/denom
    return mu


def liquid_mixture_viscosity(target):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    mus = [c['pore.viscosity'] for c in target.components.values()]
    mu = 0.0
    for i in range(len(xs)):
        mu += xs[i]*np.log(mus[i])
    return np.exp(mu)


def liquid_mixture_critical_volume(target):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    Vcs = [c.parameters['critical_volume'] for c in target.components.values()]
    N = len(xs)  # Number of components

    Vm1 = np.sum([xs[i]*Vcs[i] for i in range(N)], axis=0)
    Vm2 = np.sum([xs[i]*(Vcs[i])**(2/3) for i in range(N)], axis=0)
    Vm3 = np.sum([xs[i]*(Vcs[i])**(1/3) for i in range(N)], axis=0)
    Vm = 0.25*(Vm1 + 3*(Vm2*Vm3))
    return Vm


def liquid_mixture_critical_temperature(target, Vc='pore.critical_volume'):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    Tcs = [c.parameters['critical_temperature'] for c in target.components.values()]
    Vcs = [c.parameters['critical_volume'] for c in target.components.values()]
    Vm = target[Vc]
    N = len(xs)  # Number of components

    num = np.zeros_like(xs[0])
    for i in range(N):
        for j in range(N):
            VT = (Vcs[i]*Tcs[i]*Vcs[j]*Tcs[j])**0.5
            num += xs[i]*xs[j]*VT
    Tcm = num/Vm
    return Tcm


def mixture_acentric_factor(target):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    omegas = [c.parameters['acentric_factor'] for c in target.components.values()]
    omega = np.sum([omegas[i]*xs[i] for i in range(len(xs))], axis=0)
    return omega


def mixture_molecular_weight(target):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    MWs = [c.parameters['molecular_weight'] for c in target.components.values()]
    MW = np.zeros_like(xs[0])
    for i in range(len(xs)):
        MW += xs[i]*MWs[i]
    return MW


def liquid_mixture_density(target, temperature='pore.temperature'):
    # Actually pure component density using mixture properties
    T = target[temperature]
    MW = target['pore.molecular_weight']
    Tc = target['pore.critical_temperature']
    Vc = target['pore.critical_volume']
    omega = target['pore.acentric_factor']
    T = np.clip(T, -np.inf, Tc)
    Tr = T/Tc
    tau = 1.0 - Tr
    tau_cbrt = (tau)**(1.0/3.)
    a = 0.296123
    b = 0.0480645
    c = 0.0427258
    d = 0.386914
    e = 0.190454
    f = 0.81446
    g = 1.43907
    h = 1.52816
    V_delta = (-a + Tr*(Tr*(-b*Tr - c) + d))/(Tr - 1.00001)
    V_0 = tau_cbrt*(tau_cbrt*(tau_cbrt*(e*tau_cbrt - f) + g) - h) + 1.0
    Vm = Vc*V_0*(1.0 - omega*V_delta)
    rhoL = Vm_to_rho(Vm=Vm, MW=MW)
    return rhoL


def liquid_mixture_thermal_conductivity(target, density='pore.density'):
    xs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    kLs = [c['pore.thermal_conductivity'] for c in target.components.values()]
    # kL = numba_vectorized.DIPPR9I(xs, kLs)  # Another one that doesn't work
    Vm = [rho_to_Vm(c[density], c.parameters['molecular_weight'])
          for c in target.components.values()]
    denom = np.sum([xs[i]*Vm[i] for i in range(len(xs))], axis=0)
    phis = np.array([xs[i]*Vm[i] for i in range(len(xs))])/denom
    kij = 2/np.sum([1/kLs[i] for i in range(len(xs))], axis=0)
    kmix = np.zeros_like(xs[0])
    N = len(xs)
    for i in range(N):
        for j in range(N):
            kmix += phis[i]*phis[j]*kij
    return kmix


def gas_mixture_thermal_conductivity(target, temperature='pore.temperature'):
    T = target[temperature]
    ys = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    kGs = [c['pore.thermal_conductivity'] for c in target.components.values()]
    MWs = [c.parameters['molecular_weight'] for c in target.components.values()]
    # numba_vectorized.Lindsay_Bromley(T, ys, kGs, mus, Tbs, MWs)
    kmix = np.zeros_like(T)
    for i in range(len(ys)):
        num = ys[i]*kGs[i]
        A = 0.0
        denom = 0.0
        for j in range(len(ys)):
            A += np.sqrt(MWs[j]/MWs[i])
            denom += ys[j]*A
        kmix += num/denom
    return kmix


def mixture_heat_capacity(target):
    zs = [target['pore.mole_fraction.' + c.name] for c in target.components.values()]
    Cps = [c['pore.heat_capacity'] for c in target.components.values()]
    Cpmix = np.sum([zs[i]*Cps[i] for i in range(len(zs))], axis=0)
    return Cpmix


def gas_mixture_LJ_collision_integral(
        target,
        temperature='pore.temperature',
        epsilon='pore.LJ_epsilon',
        ):
    T = target[temperature]
    eAB = target[epsilon]
    A = 1.06036
    B = 0.15610
    C = 0.19300
    D = 0.47635
    E = 1.03587
    F = 1.52996
    G = 1.76474
    H = 3.89411
    Ts = k*T/eAB
    Omega = A/Ts**B + C/np.exp(D*Ts) + E/np.exp(F*Ts) + G/np.exp(H*Ts)
    return Omega


def gas_mixture_LJ_epsilon(target):
    es = [c.parameters['lennard_jones_epsilon'] for c in target.components.values()]
    eAB = np.sqrt(np.prod(es))
    return eAB


def gas_mixture_LJ_sigma(target):
    ss = [c.parameters['lennard_jones_sigma'] for c in target.components.values()]
    sAB = np.mean(ss)
    return sAB


def gas_mixture_diffusivity(
        target,
        temperature='pore.temperature',
        pressure='pore.pressure',
        sigma='pore.LJ_sigma',
        omega='pore.LJ_omega',
        ):
    MW = [c.parameters['molecular_weight'] for c in target.components.values()]
    MWAB = 2/np.sum(1/np.array(MW))
    T = target[temperature]
    P = target[pressure]
    sAB = target[sigma]
    Omega = target[omega]
    DAB = 0.00266*T*1.5/(P*MWAB**0.5*sAB**2*Omega)*101.325  # Convert to m2/s
    return DAB
