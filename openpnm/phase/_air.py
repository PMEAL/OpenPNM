from openpnm.phase import GenericPhase, GasMixture, GasByName
import openpnm.models as mods
from openpnm.utils import Docorator, Workspace


ws = Workspace()
docstr = Docorator()
__all__ = [
    'Air',
    'AirMixture',
    ]


@docstr.dedent
class Air(GenericPhase):
    r"""
    Creates a Phase object with preset models and values for air.

    Parameters
    ----------
    %(GenericPhase.parameters)s

    References
    ----------
    The correlations and constants for this class are taken from:

    ::

        E.W. Lemmon and R.T. Jacobsen, "Viscosity and Thermal Conductivity
        Equations for Nitrogen, Oxygen, Argon, and Air", Int. J. of
        Thermophysics, Vol. 25, No. 1, January 2004, pp. 21-69

    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self['pore.molecular_weight'] = 0.0291
        self['pore.critical_pressure'] = 3.786E6
        self['pore.critical_temperature'] = 132.5
        self['pore.critical_volume'] = 0.002917
        self['pore.contact_angle'] = 180.0
        self['pore.surface_tension'] = 0.072
        self.add_model(propname='pore.molar_density',
                       model=mods.phase.molar_density.ideal_gas)
        self.add_model(propname='pore.diffusivity',
                       model=mods.phase.diffusivity.fuller,
                       MA=0.032, MB=0.028,
                       vA=16.6, vB=17.9)
        self.add_model(propname='pore.thermal_conductivity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[0.00422791, 0.0000789606, -1.56383E-08])
        self.add_model(propname='pore.electrical_conductivity',
                       model=mods.misc.constant,
                       value=1e-15)
        self.add_model(propname='pore.viscosity',
                       model=mods.misc.polynomial,
                       prop='pore.temperature',
                       a=[0.00000182082, 6.51815E-08, -3.48553E-11,
                          1.11409E-14])


class AirMixture(GasMixture):
    r"""
    A Mixture class consisting of two species classes for O2 and N2

    """

    def __init__(self, **kwargs):
        network = kwargs['network']
        o2 = GasByName(network=network, species='O2', settings={'prefix': 'O2'})
        n2 = GasByName(network=network, species='N2', settings={'prefix': 'N2'})
        super().__init__(components=[o2, n2], **kwargs)
        self['pore.mole_fraction.' + o2.name] = 0.21
        self['pore.mole_fraction.' + n2.name] = 0.79

        self.add_model(propname='pore.molecular_weight',
                       model=mods.phases.mixtures.mole_weighted_average,
                       prop='param.molecular_weight')
        self.add_model(propname='pore.density',
                       model=mods.phases.density.ideal_gas,
                       temperature='pore.temperature',
                       pressure='pore.pressure',
                       mol_weight='pore.molecular_weight')









