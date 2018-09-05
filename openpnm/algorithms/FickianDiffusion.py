from openpnm.algorithms import ReactiveTransport
from openpnm.utils import logging, Docorator
logger = logging.getLogger(__name__)
docstr = Docorator()


@docstr.dedent
class FickianDiffusion(ReactiveTransport):
    r"""
    A class to simulate binary diffusion.

    Parameters
    ----------
    %(GenericTransport.class.parameters)s

    Notes
    -----
    Fickian diffusion in porous materials occurs in the void space, but
    because the diffusion is confined to pores it is impacted by the porosity
    and tortuosity of the network.  Thus the total diffusive flux through the
    network is reduced.  This class can be used to simualte diffusion-reaction
    in domains with customized boundary conditions, or it can be used
    to calculate the effective diffusivity of the network by applying
    controlled boundary conditions on opposing faces, calculate the diffusion
    rate, and inverting Fick's first law:

    .. math::

        D_{eff} = N_{A}L/(A\Delta C_{A})

    This class includes a method for calculating Deff automatically assuming
    appropriate boundary conditions were applied
    (``calc_effective_diffusivity``). The length and area of the domain should
    be supplied, but if they are not an attempt is made to calculate them.

    ----

    %(GenericTransport.class.notes)s

    """

    def __init__(self, settings={}, phase=None, **kwargs):
        def_set = {'phase': None,
                   'quantity': 'pore.concentration',
                   'conductance': 'throat.diffusive_conductance',
                   'gui': {'setup':        {'phase': None,
                                            'quantity': '',
                                            'conductance': ''},
                           'set_rate_BC':  {'pores': None,
                                            'values': None},
                           'set_value_BC': {'pores': None,
                                            'values': None},
                           'set_source':   {'pores': None,
                                            'propname': ''}
                           }
                   }
        super().__init__(**kwargs)
        self.settings.update(def_set)
        self.settings.update(settings)
        if phase is not None:
            self.setup(phase=phase)

    @docstr.dedent
    def setup(self, phase=None, quantity='', conductance='', **kwargs):
        r"""
        This method takes several arguments that are essential to running the
        algorithm and adds them to the settings.

        Parameters
        ----------
        %(ReactiveTransport.setup.parameters)s

        ----
        The following settings are used by the source term iterations:

        %(ReactiveTransport.setup.other_parameters)s

        %(ReactiveTransport.setup.notes)s

        ----

        The following settings are used to control the behavior of the solver:

        %(GenericTransport.setup.other_parameters)s

        """
        if phase:
            self.settings['phase'] = phase.name
        if quantity:
            self.settings['quantity'] = quantity
        if conductance:
            self.settings['conductance'] = conductance
        super().setup(**kwargs)

    @docstr.dedent
    def calc_effective_diffusivity(self, inlets=None, outlets=None,
                                   domain_area=None, domain_length=None):
        r"""
        Calculates the effective diffusivity of the network

        Parameters
        ----------
        %(GenericTransport._calc_eff_prop.parameters)s

        Notes
        -----
        %(GenericTransport._calc_eff_prop.notes)s

        """
        return self._calc_eff_prop(inlets=inlets, outlets=outlets,
                                   domain_area=domain_area,
                                   domain_length=domain_length)
