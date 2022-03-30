import logging
from openpnm.utils import Workspace, SettingsAttr, Docorator
from openpnm.phase import GenericPhase
from openpnm.physics import GenericPhysics
from openpnm.algorithms import FickianDiffusion
from openpnm.metrics import GenericTransportMetrics
from openpnm import models
logger = logging.getLogger(__name__)
ws = Workspace()
docstr = Docorator()


@docstr.get_sections(base='EffectiveDiffusivitySettings',
                     sections=['Parameters'])
@docstr.dedent
class EffectiveDiffusivitySettings:
    r"""
    Defines the settings for EffectiveDiffusivity

    ----------
    prefix : str
        The default prefix to use when generating a name
    inlet : str
        The pore labels for diffusion inlet.
    outlet : str
        The pore labels for diffusion outlet.
    area : scalar
        The cross sectional area of the network relative to the inlet and outlet
    length: scalar
        The length of the network relative to the inlet and outlet

    """
    prefix = 'edif'
    inlet = 'left'
    outlet = 'right'
    area = None
    length = None


class EffectiveDiffusivity(GenericTransportMetrics):
    r"""
    This class works by applying 'value' boundary conditions across the
    domain to find molar flow, then using Fick's law to back-calculate
    the effective diffusivity of the domain.

    .. math::

        D_{eff} = \frac{n_{A} L }{\Delta C_{A} A }

    and

    .. math::

        D_{eff} = D_{AB} \frac{\varepsilon}{\tau}


    Examples
    --------
    >>> import openpnm as op
    >>> import numpy as np
    >>> np.random.seed(5)
    >>> pn = op.network.Cubic(shape=[10, 10, 10], spacing=1e-5)
    >>> geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)

    Now find the effective diffusivity of the network:

    >>> ED = op.metrics.EffectiveDiffusivity(network=pn)
    >>> Deff = ED.run()
    >>> print(np.round(Deff,decimals=3))
    0.049

    """

    def __init__(self, settings=None, **kwargs):
        self.settings = SettingsAttr(EffectiveDiffusivitySettings, settings)
        super().__init__(settings=self.settings, **kwargs)

    def run(self, verbose=False):
        r"""
        Execute the diffusion simulations in the principle directions.

        """
        phase = GenericPhase(network=self.network)
        phase['pore.diffusivity'] = 1.0
        phase['throat.diffusivity'] = 1.0
        mod = models.physics.diffusive_conductance.ordinary_diffusion
        for geom in self.project.geometries().values():
            phys = GenericPhysics(network=self.network,
                                  phase=phase, geometry=geom)
            phys.add_model(propname='throat.diffusive_conductance', model=mod)
        inlet = self.network.pores(self.settings['inlet'])
        outlet = self.network.pores(self.settings['outlet'])
        Diff = FickianDiffusion(network=self.project.network, phase=phase)
        Diff.set_value_BC(pores=inlet, values=1.0)
        Diff.set_value_BC(pores=outlet, values=0.0)
        Diff.run(verbose=verbose)
        Deff = self._calc_eff_prop(inlets=inlet, outlets=outlet,
                                   domain_area=self.settings['area'],
                                   domain_length=self.settings['length'],
                                   rates=Diff.rate(pores=inlet),
                                   prop_diff=1)
        # Deff = R*L/A  # Conc gradient and diffusivity were both unity
        return Deff
