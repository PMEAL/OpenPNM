from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator, SettingsAttr
docstr = Docorator()
logger = logging.getLogger(__name__)


@docstr.get_sections(base='FickianDiffusionSettings',
                     sections=['Parameters'])
@docstr.dedent
class FickianDiffusionSettings():
    r"""

    Parameters
    ----------
    %(GenericTransportSettings.parameters)s
    quantity : str (default = 'pore.concentration')
        The name of the physical quantity to be calculated
    conductance : str (default = 'throat.diffusive_conductance')
        The name of the pore-scale transport conductance values. These are
        typically calculated by a model attached to a *Physics* object
        associated with the given *Phase*.

    Other Parameters
    ----------------

    **The following parameters pertain to the ReactiveTransport class**

    %(ReactiveTransportSettings.other_parameters)s

    ----

    **The following parameters pertain to the GenericTransport class**

    %(GenericTransportSettings.other_parameters)s

    """
    prefix = 'fick'
    quantity = 'pore.concentration'
    conductance = 'throat.diffusive_conductance'


@docstr.dedent
class FickianDiffusion(ReactiveTransport):
    r"""
    A class to simulate binary diffusion with reactions

    Parameters
    ----------
    %(ReactiveTransport.parameters)s

    Notes
    -----
    Fickian diffusion in porous materials occurs in the void space, but
    becuase the diffusion is defined to pores it is impacted by the porosity
    and tortuosity of the network.  Thus the total diffusive flux through the
    network is reduced.  This class can be used to simualte diffusion-reaction
    in domains with arbitrarily complex boundary conditions, or it can be used
    to calculate the effective diffusivity of the network by applying
    controlled boundary conditions on opposing faces, calculate the diffusion
    rate, and inverting Fick's first law:

    .. math::

        D_{eff} = N_{A}*L/(A*\Delta C_{A})

    This class includes a method for calculating Deff automatically assuming
    appropriate boundary conditions were applied (``calc_eff_diffusivity``).
    The length and area of the domain should be supplied, but if they are
    not an attempt is made to calculate them.

    """

    def __init__(self, settings={}, **kwargs):
        self.settings = SettingsAttr(FickianDiffusionSettings, settings)
        super().__init__(settings=self.settings, **kwargs)

    def calc_effective_diffusivity(self, inlets=None, outlets=None,
                                   domain_area=None, domain_length=None):
        r"""
        This calculates the effective diffusivity in this linear transport
        algorithm.

        Parameters
        ----------
        inlets : array_like
            The pores where the inlet composition boundary conditions were
            applied.  If not given an attempt is made to infer them from the
            algorithm.

        outlets : array_like
            The pores where the outlet composition boundary conditions were
            applied.  If not given an attempt is made to infer them from the
            algorithm.

        domain_area : scalar, optional
            The area of the inlet (and outlet) boundary faces.  If not given
            then an attempt is made to estimate it, but it is usually
            underestimated.

        domain_length : scalar, optional
            The length of the domain between the inlet and outlet boundary
            faces.  If not given then an attempt is made to estimate it, but it
            is usually underestimated.

        Notes
        -----
        The area and length of the domain are found using the bounding box
        around the inlet and outlet pores which do not necessarily lie on the
        edge of the domain, resulting in underestimation of sizes.
        """
        return self._calc_eff_prop(inlets=inlets, outlets=outlets,
                                   domain_area=domain_area,
                                   domain_length=domain_length)


if __name__ == "__main__":
    import openpnm as op
    pn = op.network.Cubic([3, 3, 3])
    geo = op.geometry.SpheresAndCylinders(network=pn, pores=pn.Ps, throats=pn.Ts)
    air = op.phases.Air(network=pn)
    phys = op.physics.Basic(network=pn, phase=air, geometry=geo)
    s = {'bob': 3}
    fd = op.algorithms.FickianDiffusion(network=pn, phase=air, settings=s)
    print(fd.settings)
